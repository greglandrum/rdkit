//
//  Copyright (C) 2004-2019 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY

#include <RDBoost/python.h>

#include <string>
#include <memory>
#ifdef RDK_BUILD_THREADSAFE_SSS

#include <GraphMol/FileParsers/GeneralFileReader.h>
#include <GraphMol/RDKitBase.h>

namespace python = boost::python;

namespace RDKit {
namespace {
python::tuple determineFormatHelper(const std::string &fileName){
  std::string fmt="";
  std::string cfmt="";
  GeneralMolSupplier::determineFormat(fileName,fmt,cfmt);
  return python::make_tuple(fmt,cfmt);
}

}  // namespace


struct generalizedfilereader_wrap {
  static void wrap() {
    python::def("DetermineFileFormat",determineFormatHelper,"returns a two-tuple with (file_format, compression_format)");
  }
};
}  // namespace RDKit

void wrap_generalizedFileReader() { RDKit::generalizedfilereader_wrap::wrap(); }
#endif