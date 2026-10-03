//
//  Copyright (C) 2025 and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/AvalonGenerator.h>

using namespace RDKit;
namespace nb = nanobind;
using namespace nb::literals;

namespace RDKit {
namespace AvalonWrapper {

void exportAvalon(nb::module_ &m) {
  nb::class_<AvalonFP::AvalonArguments, FingerprintArguments>(
      m, "AvalonFingerprintOptions")
      .def_rw("bitFlags", &AvalonFP::AvalonArguments::d_bitFlags,
              "Avalon feature flags (USE_* in the Avalon toolkit)")
      .def_rw("isQuery", &AvalonFP::AvalonArguments::df_isQuery,
              "generate the fingerprint for a query molecule")
      .def_rw("accumulateAsQuery",
              &AvalonFP::AvalonArguments::df_accumulateAsQuery,
              "use query mode in the second (Daylight aromaticity) pass; "
              "matches AvalonTools.GetAvalonCountFP()");

  m.def(
      "GetAvalonGenerator",
      [](std::uint32_t fpSize, std::uint32_t bitFlags, bool isQuery,
         bool accumulateAsQuery) {
        return AvalonFP::getAvalonGenerator<std::uint32_t>(
            fpSize, bitFlags, isQuery, accumulateAsQuery);
      },
      "fpSize"_a = 512, "bitFlags"_a = AvalonFP::defaultAvalonBitFlags,
      "isQuery"_a = false, "accumulateAsQuery"_a = false,
      R"DOC(Get an Avalon fingerprint generator

ARGUMENTS:
    - fpSize: size of the generated fingerprint
    - bitFlags: Avalon feature flags
    - isQuery: generate the fingerprint for a query molecule
    - accumulateAsQuery: use query mode in the second pass (set this to
      reproduce AvalonTools.GetAvalonCountFP())

Bit vectors correspond to Avalon's SetFingerprintBits(), count fingerprints
to SetFingerprintCountsWithFocus().
fromAtoms and ignoreAtoms are not supported.

RETURNS: FingerprintGenerator
)DOC",
      nb::rv_policy::take_ownership);
}
}  // namespace AvalonWrapper
}  // namespace RDKit
