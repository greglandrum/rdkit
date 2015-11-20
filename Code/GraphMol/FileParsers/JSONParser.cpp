//
//  Copyright (C) 2015 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "FileParsers.h"
#include "FileParserUtils.h"
#include "MolFileStereochem.h"
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/RDLog.h>

#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/LocaleSwitcher.h>
#include <typeinfo>
#include <exception>
#include <sstream>
#include <locale>
#include <stdlib.h>
#include <rapidjson/document.h>
#include <rapidjson/error/en.h>

namespace RDKit {

namespace JSONParserUtils {

// from the rapidjson documentation
class IStreamWrapper {
 public:
  typedef char Ch;
  IStreamWrapper(std::istream &is) : is_(is) {}
  Ch Peek() const {  // 1
    int c = is_.peek();
    return c == std::char_traits<char>::eof() ? '\0' : (Ch)c;
  }
  Ch Take() {  // 2
    int c = is_.get();
    return c == std::char_traits<char>::eof() ? '\0' : (Ch)c;
  }
  size_t Tell() const { return (size_t)is_.tellg(); }  // 3
  Ch *PutBegin() {
    assert(false);
    return 0;
  }
  void Put(Ch) { assert(false); }
  void Flush() { assert(false); }
  size_t PutEnd(Ch *) {
    assert(false);
    return 0;
  }

 private:
  IStreamWrapper(const IStreamWrapper &);
  IStreamWrapper &operator=(const IStreamWrapper &);
  std::istream &is_;
};

void ParseJSONAtoms(std::istream *inStream, unsigned int &line,
                    unsigned int nAtoms, RWMol *mol, Conformer *conf) {
  PRECONDITION(inStream, "bad stream");
  PRECONDITION(mol, "bad molecule");
  PRECONDITION(conf, "bad conformer");
  for (unsigned int i = 0; i < nAtoms; ++i) {
    RDGeom::Point3D pos;
    Atom *atom = NULL;
    unsigned int aid = mol->addAtom(atom, false, true);
    conf->setAtomPos(aid, pos);
  }
}

// returns whether or not any sign of chirality was detected
void ParseJSONBonds(std::istream *inStream, unsigned int &line,
                    unsigned int nBonds, RWMol *mol, bool &chiralityPossible) {
  PRECONDITION(inStream, "bad stream");
  PRECONDITION(mol, "bad molecule");
  for (unsigned int i = 0; i < nBonds; ++i) {
    Bond *bond = NULL;
    // if we got an aromatic bond set the flag on the bond and the connected
    // atoms
    if (bond->getBondType() == Bond::AROMATIC) {
      bond->setIsAromatic(true);
      mol->getAtomWithIdx(bond->getBeginAtomIdx())->setIsAromatic(true);
      mol->getAtomWithIdx(bond->getEndAtomIdx())->setIsAromatic(true);
    }
    // if the bond might have chirality info associated with it, set a flag:
    if (bond->getBondDir() != Bond::NONE &&
        bond->getBondDir() != Bond::UNKNOWN) {
      chiralityPossible = true;
    }
    mol->addBond(bond, true);
  }
}

bool ParseJSONProperties(std::istream *inStream, unsigned int &line,
                         RWMol *mol) {
  PRECONDITION(inStream, "bad stream");
  PRECONDITION(mol, "bad molecule");

  bool fileComplete = false;
  return fileComplete;
}

namespace {
int getIntDefaultValue(const char *key, const rapidjson::Value &from,
                       const rapidjson::Value &defaults) {
  rapidjson::Value::ConstMemberIterator miter = from.FindMember(key);
  if (miter != from.MemberEnd()) {
    if (!miter->value.IsInt())
      throw FileParseException(std::string("Bad format: ") + std::string(key) +
                               std::string(" is not an int"));
    return miter->value.GetInt();
  } else {
    rapidjson::Value::ConstMemberIterator miter = defaults.FindMember(key);
    if (miter != defaults.MemberEnd()) {
      if (!miter->value.IsInt())
        throw FileParseException(std::string("Bad format: default value of ") +
                                 std::string(key) +
                                 std::string(" is not an int"));
      return miter->value.GetInt();
    }
  }
  return 0;
}
}  // end of anonymous namespace

RWMol *JSONDocumentToMol(rapidjson::Document &jsondoc, bool sanitize,
                         bool removeHs, bool strictParsing) {
  // TODO:
  //  - stereochemistry
  //  - representations
  //
  std::string tempStr;
  bool fileComplete = false;
  bool chiralityPossible = false;
  Utils::LocaleSwitcher ls;

  // some error checking
  if (!jsondoc.IsObject())
    throw FileParseException("Bad Format: JSON should be an object");
  if (!jsondoc.HasMember("type") || jsondoc["type"] != "mol")
    throw FileParseException("Bad Format: missing or bad type in JSON");

  RWMol *res = new RWMol();
  Conformer *conf = NULL;
  try {
    if (jsondoc.HasMember("title")) {
      if (!jsondoc["title"].IsString())
        throw FileParseException("Bad Format: title should be a string");
      res->setProp(common_properties::_Name, jsondoc["title"].GetString());
    }
  } catch (FileParseException &e) {
    // catch our exceptions and throw them back after cleanup
    delete res;
    delete conf;
    res = NULL;
    conf = NULL;
    throw e;
  }

  fileComplete = true;

  int dimension = 2;
  if (jsondoc.HasMember("dimension")) {
    if (!jsondoc["dimension"].IsNumber())
      throw FileParseException("Bad Format: dimension should be int");
    dimension = jsondoc["dimension"].GetInt();
  }

  // -----------
  // parse atoms:
  rapidjson::Value atomDefaults;
  if (jsondoc.HasMember("atomdefaults")) {
    atomDefaults = jsondoc["atomdefaults"];
    if (!atomDefaults.IsObject())
      throw FileParseException("Bad Format: atomdefaults is not an object");
  }

  size_t nAtoms = 0;
  if (jsondoc.HasMember("atoms")) {
    const rapidjson::Value &atoms = jsondoc["atoms"];
    if (!atoms.IsArray())
      throw FileParseException("Bad Format: atoms not in a json array");
    nAtoms = atoms.Size();
    conf = new Conformer(nAtoms);
    if (dimension == 3) conf->set3D(true);
    for (size_t i = 0; i < nAtoms; ++i) {
      const rapidjson::Value &av = atoms[i];
      unsigned int num = getIntDefaultValue("element", av, atomDefaults);
      Atom *atom = new Atom(num);

      atom->setNumExplicitHs(
          getIntDefaultValue("implicithcount", av, atomDefaults));
      atom->setNoImplicit(true);
      atom->setFormalCharge(
          getIntDefaultValue("formalcharge", av, atomDefaults));

      if (av.HasMember("coords")) {
        const rapidjson::Value &cv = av["coords"];
        if (!cv.IsArray())
          throw FileParseException("Bad Format: coords not in a json array");
        if (cv.Size() < dimension)
          throw FileParseException("Bad Format: coords array too short");
        RDGeom::Point3D p(0, 0, 0);
        p.x = cv[0].GetDouble();
        p.y = cv[1].GetDouble();
        if (dimension > 2) p.z = cv[2].GetDouble();
        conf->setAtomPos(i, p);
      }

      // add the atom
      res->addAtom(atom, true, true);
    }
    res->addConformer(conf);
    conf = NULL;
  }

  // -----------
  // parse bonds:
  rapidjson::Value bondDefs;
  if (jsondoc.HasMember("bonddefaults")) {
    bondDefs = jsondoc["bonddefaults"];
    if (!bondDefs.IsObject())
      throw FileParseException("Bad Format: bonddefaults is not an object");
  }

  size_t nBonds = 0;
  if (jsondoc.HasMember("bonds")) {
    const rapidjson::Value &bonds = jsondoc["bonds"];
    if (!bonds.IsArray())
      throw FileParseException("Bad Format: bonds not in a json array");
    nBonds = bonds.Size();

    for (size_t i = 0; i < nBonds; ++i) {
      const rapidjson::Value &bv = bonds[i];
      unsigned int oi = getIntDefaultValue("order", bv, bondDefs);
      Bond::BondType bt;
      switch (oi) {
        case 0:
          bt = Bond::UNSPECIFIED;
          break;
        case 1:
          bt = Bond::SINGLE;
          break;
        case 2:
          bt = Bond::DOUBLE;
          break;
        case 3:
          bt = Bond::TRIPLE;
          break;
        default:
          throw FileParseException("Bad Format: unrecognized bond order");
      }

      if (!bv.HasMember("atoms"))
        throw FileParseException(
            "Bad Format: bond object has no atoms element");
      const rapidjson::Value &bAts = bv["atoms"];
      if (!bonds.IsArray())
        throw FileParseException("Bad Format: bond atoms not in a json array");
      if (bAts.Size() != 2)
        throw FileParseException("Bad Format: bond atoms array not 2 long");
      unsigned int begAt = bAts[0].GetInt();
      unsigned int endAt = bAts[1].GetInt();

      // add the bond
      unsigned int bidx = res->addBond(begAt, endAt, bt);
      // do some stuff
    }
  }

  if (!fileComplete) {
    delete res;
    delete conf;
    res = NULL;
    conf = NULL;
    throw FileParseException("file incomplete");
  }

  if (res) {
    // calculate explicit valence on each atom:
    for (RWMol::AtomIterator atomIt = res->beginAtoms();
         atomIt != res->endAtoms(); ++atomIt) {
      (*atomIt)->calcExplicitValence(false);
    }

    // postprocess mol file flags
    // ProcessMolProps(res);

    // update the chirality and stereo-chemistry
    //
    // NOTE: we detect the stereochemistry before sanitizing/removing
    // hydrogens because the removal of H atoms may actually remove
    // the wedged bond from the molecule.  This wipes out the only
    // sign that chirality ever existed and makes us sad... so first
    // perceive chirality, then remove the Hs and sanitize.
    //
    // One exception to this (of course, there's always an exception):
    // DetectAtomStereoChemistry() needs to check the number of
    // implicit hydrogens on atoms to detect if things can be
    // chiral. However, if we ask for the number of implicit Hs before
    // we've called MolOps::cleanUp() on the molecule, we'll get
    // exceptions for common "weird" cases like a nitro group
    // mis-represented as -N(=O)=O.  *SO*... we need to call
    // cleanUp(), then detect the stereochemistry.
    // (this was Issue 148)
    //
    const Conformer &conf = res->getConformer();
    if (chiralityPossible) {
      MolOps::cleanUp(*res);
      DetectAtomStereoChemistry(*res, &conf);
    }

    if (sanitize) {
      try {
        if (removeHs) {
          MolOps::removeHs(*res, false, false);
        } else {
          MolOps::sanitizeMol(*res);
        }

        // now that atom stereochem has been perceived, the wedging
        // information is no longer needed, so we clear
        // single bond dir flags:
        ClearSingleBondDirFlags(*res);

        // unlike DetectAtomStereoChemistry we call DetectBondStereoChemistry
        // here after sanitization because we need the ring information:
        DetectBondStereoChemistry(*res, &conf);
      } catch (...) {
        delete res;
        res = NULL;
        throw;
      }
      MolOps::assignStereochemistry(*res, true, true, true);
    } else {
      // we still need to do something about double bond stereochemistry
      // (was github issue 337)
      DetectBondStereoChemistry(*res, &conf);
    }
  }
  return res;
}

}  // end of JSONParserUtils namespace

//------------------------------------------------
//
//  Read a molecule from a stream
//
//------------------------------------------------
RWMol *JSONDataStreamToMol(std::istream *inStream, bool sanitize, bool removeHs,
                           bool strictParsing) {
  PRECONDITION(inStream, "no stream");

  rapidjson::Document jsondoc;
  JSONParserUtils::IStreamWrapper isw(*inStream);
  if (jsondoc.ParseStream(isw).HasParseError()) {
    std::string msg = rapidjson::GetParseError_En(jsondoc.GetParseError());
    BOOST_LOG(rdErrorLog) << "JSON Parse Error: " << msg
                          << " at position: " << jsondoc.GetErrorOffset()
                          << std::endl;
    throw FileParseException(msg);
  }
  return JSONParserUtils::JSONDocumentToMol(jsondoc, sanitize, removeHs,
                                            strictParsing);
};

RWMol *JSONDataStreamToMol(std::istream &inStream, bool sanitize, bool removeHs,
                           bool strictParsing) {
  return JSONDataStreamToMol(&inStream, sanitize, removeHs, strictParsing);
};
//------------------------------------------------
//
//  Read a molecule from a string
//
//------------------------------------------------
RWMol *JSONToMol(const std::string &json, bool sanitize, bool removeHs,
                 bool strictParsing) {
  rapidjson::Document jsondoc;
  if (jsondoc.Parse(json.c_str()).HasParseError()) {
    std::string msg = rapidjson::GetParseError_En(jsondoc.GetParseError());
    BOOST_LOG(rdErrorLog) << "JSON Parse Error: " << msg
                          << " at position: " << jsondoc.GetErrorOffset()
                          << std::endl;
    throw FileParseException(msg);
  }
  return JSONParserUtils::JSONDocumentToMol(jsondoc, sanitize, removeHs,
                                            strictParsing);
}
}
