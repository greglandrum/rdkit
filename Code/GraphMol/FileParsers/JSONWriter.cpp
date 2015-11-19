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

std::string MolToJSON(const ROMol &mol, bool includeStereo, int confId,
                      bool kekulize) {
  rj::Document d;

  rj::StringBuffer buffer;
  rj::Writer<rj::StringBuffer> writer(buffer);
  d.Accept(writer);
  return std::string(buffer.GetString());
}
}
