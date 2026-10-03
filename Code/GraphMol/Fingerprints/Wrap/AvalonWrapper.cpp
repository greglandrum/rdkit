//
//  Copyright (C) 2025 and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <boost/python.hpp>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/AvalonGenerator.h>
#include <RDBoost/Wrap.h>

using namespace RDKit;
namespace python = boost::python;

namespace RDKit {
namespace AvalonWrapper {

FingerprintGenerator<std::uint32_t> *getAvalonGen(std::uint32_t fpSize,
                                                  std::uint32_t bitFlags,
                                                  bool isQuery,
                                                  bool accumulateAsQuery) {
  return AvalonFP::getAvalonGenerator<std::uint32_t>(fpSize, bitFlags, isQuery,
                                                     accumulateAsQuery);
}

void exportAvalon() {
  python::class_<AvalonFP::AvalonArguments, python::bases<FingerprintArguments>,
                 boost::noncopyable>("AvalonFingerprintOptions",
                                     python::no_init)
      .def_readwrite("bitFlags", &AvalonFP::AvalonArguments::d_bitFlags,
                     "Avalon feature flags (USE_* in the Avalon toolkit)")
      .def_readwrite("isQuery", &AvalonFP::AvalonArguments::df_isQuery,
                     "generate the fingerprint for a query molecule")
      .def_readwrite("accumulateAsQuery",
                     &AvalonFP::AvalonArguments::df_accumulateAsQuery,
                     "use query mode in the second (Daylight aromaticity) "
                     "pass; matches AvalonTools.GetAvalonCountFP()");
  python::def(
      "GetAvalonGenerator", &getAvalonGen,
      (python::arg("fpSize") = 512,
       python::arg("bitFlags") = AvalonFP::defaultAvalonBitFlags,
       python::arg("isQuery") = false,
       python::arg("accumulateAsQuery") = false),
      "Get an Avalon fingerprint generator\n\n"
      "  ARGUMENTS:\n"
      "    - fpSize: size of the generated fingerprint\n"
      "    - bitFlags: Avalon feature flags\n"
      "    - isQuery: generate the fingerprint for a query molecule\n"
      "    - accumulateAsQuery: use query mode in the second pass (set this "
      "to reproduce AvalonTools.GetAvalonCountFP())\n\n"
      "Bit vectors correspond to Avalon's SetFingerprintBits(), count "
      "fingerprints to SetFingerprintCountsWithFocus().\n"
      "fromAtoms and ignoreAtoms are not supported.\n\n"
      "  RETURNS: FingerprintGenerator\n\n",
      python::return_value_policy<python::manage_new_object>());
}
}  // namespace AvalonWrapper
}  // namespace RDKit
