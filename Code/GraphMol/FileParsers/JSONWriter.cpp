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

#include <fstream>
#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

namespace rj = rapidjson;

namespace RDKit {

namespace JSONWriterUtils {}  // end of JSONParserUtils namespace

//------------------------------------------------
//
//  Generate JSON from a molecule
//
//------------------------------------------------

std::string MolToJSON(const ROMol& mol, bool includeStereo, int confId,
                      bool kekulize) {
  const Conformer* conf = NULL;
  if (mol.getNumConformers()) conf = &mol.getConformer(confId);

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
  mol.getPropIfPresent(common_properties::_Name, nm);
  val = rj::StringRef(nm.c_str());
  doc.AddMember("title", val, allocator);

  // -----
  // write atoms
  // FIX: defaults
  rj::Value atoms(rj::kArrayType);
  atoms.Reserve(mol.getNumAtoms(), allocator);
  for (RWMol::ConstAtomIterator at = mol.beginAtoms(); at != mol.endAtoms();
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
  bonds.Reserve(mol.getNumBonds(), allocator);
  for (RWMol::ConstBondIterator bnd = mol.beginBonds(); bnd != mol.endBonds();
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

    bonds.PushBack(bndV, allocator);
  }
  doc.AddMember("bonds", bonds, allocator);

  // get the text and return it
  rj::StringBuffer buffer;
  rj::Writer<rj::StringBuffer> writer(buffer);
  doc.Accept(writer);
  // FIX: ensure this doesn't leak
  return std::string(buffer.GetString());
}
}
