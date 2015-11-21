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
#include "MolFileStereochem.h"
#include <RDGeneral/RDLog.h>
#include <RDGeneral/types.h>
#include <RDGeneral/versions.h>

#include <fstream>
#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

namespace rj = rapidjson;
namespace RDKit {
namespace JSONWriterUtils {}  // end of JSONParserUtils namespace

//------------------------------------------------
//
//  Generate JSON from a molecule.
//
//------------------------------------------------

std::string MolToJSON(const ROMol& mol, int confId, bool kekulize) {
  ROMol tromol(mol);
  RWMol& trwmol = static_cast<RWMol&>(tromol);
  // NOTE: kekulize the molecule before writing it out
  // because of the way mol files handle aromaticity
  if (trwmol.needsUpdatePropertyCache()) {
    trwmol.updatePropertyCache(false);
  }
  if (kekulize) MolOps::Kekulize(trwmol, false);

  const Conformer* conf = NULL;
  if (trwmol.getNumConformers()) conf = &trwmol.getConformer(confId);

  rj::Document doc;
  doc.SetObject();
  rj::Document::AllocatorType& allocator = doc.GetAllocator();
  rj::Value val;
  val = "mol";
  doc.AddMember("type", val, allocator);
  val = "0.9";
  doc.AddMember("version", val, allocator);
  val = (conf && conf->is3D()) ? 3 : 2;
  doc.AddMember("dimension", val, allocator);
  std::string nm = "";
  trwmol.getPropIfPresent(common_properties::_Name, nm);
  val = rj::StringRef(nm.c_str());
  doc.AddMember("title", val, allocator);
  // -----
  // write atoms
  // FIX: defaults
  rj::Value atoms(rj::kArrayType);
  atoms.Reserve(trwmol.getNumAtoms(), allocator);
  for (RWMol::AtomIterator at = trwmol.beginAtoms(); at != trwmol.endAtoms();
       ++at) {
    rj::Value atomV(rj::kObjectType);
    val = (*at)->getAtomicNum();
    atomV.AddMember("element", val, allocator);
    rj::Value coords(rj::kArrayType);
    if (conf) {
      const RDGeom::Point3D& pos = conf->getAtomPos((*at)->getIdx());
      val = pos.x;
      coords.PushBack(val, allocator);
      val = pos.y;
      coords.PushBack(val, allocator);
      if (conf->is3D()) {
        val = pos.z;
        coords.PushBack(val, allocator);
      }
    } else {
      val = 0.0;
      coords.PushBack(val, allocator);
      val = 0.0;
      coords.PushBack(val, allocator);
    }
    atomV.AddMember("coords", coords, allocator);

    val = (*at)->getTotalNumHs(false);
    atomV.AddMember("implicithcount", val, allocator);

    val = (*at)->getFormalCharge();
    atomV.AddMember("formalcharge", val, allocator);

    val = "Undefined";  // FIX
    atomV.AddMember("stereo", val, allocator);

    atoms.PushBack(atomV, allocator);
  }
  doc.AddMember("atoms", atoms, allocator);

  // -----
  // write bonds
  // FIX: defaults
  rj::Value bonds(rj::kArrayType);
  bonds.Reserve(trwmol.getNumBonds(), allocator);
  for (RWMol::BondIterator bnd = trwmol.beginBonds(); bnd != trwmol.endBonds();
       ++bnd) {
    rj::Value bndV(rj::kObjectType);
    rj::Value atV(rj::kArrayType);
    val = (*bnd)->getBeginAtomIdx();
    atV.PushBack(val, allocator);
    val = (*bnd)->getEndAtomIdx();
    atV.PushBack(val, allocator);
    bndV.AddMember("atoms", atV, allocator);

    switch ((*bnd)->getBondType()) {
      case Bond::SINGLE:
        val = 1;
        break;
      case Bond::DOUBLE:
        val = 2;
        break;
      case Bond::TRIPLE:
        val = 3;
        break;
      default:
        val = 0;
    }
    bndV.AddMember("order", val, allocator);
    val = "Undefined";  // FIX
    bndV.AddMember("stereo", val, allocator);

    bonds.PushBack(bndV, allocator);
  }
  doc.AddMember("bonds", bonds, allocator);

  // -------
  // Now capture our internal representation
  /*
          {
              "tookit":"rdkit",
              "version":"2015.09.1",
              "aromatic_bonds":[0,1,2,3,4,5],
              "bond_rings":[ [0,1,2,3,4,5] ]
          }
  */
  rj::Value rdk(rj::kObjectType);
  val = "rdkit";
  rdk.AddMember("toolkit", val, allocator);
  val = rj::StringRef(rdkitVersion);
  rdk.AddMember("version", val, allocator);
  {
    rj::Value tarr(rj::kArrayType);
    for (unsigned int i = 0; i < trwmol.getNumBonds(); ++i) {
      if (trwmol.getBondWithIdx(i)->getIsAromatic()) {
        val = i;
        tarr.PushBack(val, allocator);
      }
    }
    rdk.AddMember("aromatic_bonds", tarr, allocator);
  }
  {
    rj::Value tarr(rj::kArrayType);
    const RingInfo* ri = trwmol.getRingInfo();
    BOOST_FOREACH (const INT_VECT& iv, ri->atomRings()) {
      rj::Value barr(rj::kArrayType);
      BOOST_FOREACH (int v, iv) {
        val = v;
        barr.PushBack(val, allocator);
      }
      tarr.PushBack(barr, allocator);
    }
    rdk.AddMember("atom_rings", tarr, allocator);
  }
  rj::Value reprs(rj::kArrayType);
  reprs.PushBack(rdk, allocator);
  doc.AddMember("representations", reprs, allocator);

  // get the text and return it
  rj::StringBuffer buffer;
  rj::Writer<rj::StringBuffer> writer(buffer);
  doc.Accept(writer);
  // FIX: ensure this doesn't leak
  return std::string(buffer.GetString());
}
}
