
#include <boost/python.hpp>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/AtomPairGenerator.h>

using namespace RDKit;
using namespace RDKit::AtomPair;
namespace python = boost::python;

namespace RDKit {
namespace AtomPairWrapper {
template <typename OutputType>
FingerprintGenerator<OutputType> *getAtomPairGenerator(
    const unsigned int minDistance, const unsigned int maxDistance,
    const bool includeChirality, const bool use2D,
    const bool useCountSimulation, python::object &py_countBounds,
    const std::uint32_t foldedSize, python::object &py_atomInvGen) {
  AtomInvariantsGenerator *atomInvariantsGenerator = nullptr;

  python::extract<AtomInvariantsGenerator *> atomInvGen(py_atomInvGen);
  if (atomInvGen.check() && atomInvGen()) {
    atomInvariantsGenerator = atomInvGen()->clone();
  }

  std::vector<std::uint32_t> countBounds = {1, 2, 4, 8};
  python::extract<std::vector<std::uint32_t>> countBoundsE(py_countBounds);
  if (countBoundsE.check() && !countBoundsE().empty()) {
    countBounds = countBoundsE();
  }

  const std::vector<std::uint32_t> countBoundsC = countBounds;

  return AtomPair::getAtomPairGenerator<OutputType>(
      minDistance, maxDistance, includeChirality, use2D,
      atomInvariantsGenerator, useCountSimulation, foldedSize, countBoundsC,
      true);
}

AtomInvariantsGenerator *getAtomPairAtomInvGen(const bool includeChirality) {
  return new AtomPairAtomInvGenerator(includeChirality);
}

void exportAtompair() {
  /*python::def(
      "GetAtomPairGenerator", &getAtomPairGenerator<std::uint32_t>,
      (python::arg("minDistance") = 1,
       python::arg("maxDistance") = AtomPair::maxPathLen - 1,
       python::arg("includeChirality") = false, python::arg("use2D") = true,
       python::arg("useCountSimulation") = true,
       python::arg("countBounds") = python::object(),
       python::arg("foldedSize") = 2048,
       python::arg("atomInvariantsGenerator") = python::object()),
      docString.c_str(),
      python::return_value_policy<python::manage_new_object>());*/

  python::def(
      "GetAtomPairGenerator", &getAtomPairGenerator<std::uint64_t>,
      (python::arg("minDistance") = 1,
       python::arg("maxDistance") = AtomPair::maxPathLen - 1,
       python::arg("includeChirality") = false, python::arg("use2D") = true,
       python::arg("useCountSimulation") = true,
       python::arg("countBounds") = python::object(),
       python::arg("foldedSize") = 2048,
       python::arg("atomInvariantsGenerator") = python::object()),
      "Get an atom pair fingerprint generator\n\n"
      "  ARGUMENTS:\n"
      "    - minDistance: minimum distance between atoms to be considered in a "
      "pair, default is 1 bond\n"
      "    - maxDistance: maximum distance between atoms to be considered in a "
      "pair, default is maxPathLen-1 bonds\n"
      "    - includeChirality: if set, chirality will be used in the atom  "
      "invariants, this is ignored if atomInvariantsGenerator is provided\n"
      "    - use2D: if set, the 2D (topological) distance matrix  will be "
      "used\n"
      "    - useCountSimulation:  if set, use count simulation while  "
      "generating the fingerprint\n"
      "    - countBounds: boundaries for count simulation, corresponding bit "
      "will be  set if the count is higher than the number provided for that "
      "spot\n"
      "    - foldedSize: size of the folded version of the fingerprints\n"
      "    - atomInvariantsGenerator: atom invariants to be used during "
      "fingerprint generation\n\n"
      "  RETURNS: FingerprintGenerator\n\n",
      python::return_value_policy<python::manage_new_object>());

  python::def("GetAtomPairAtomInvGen", &getAtomPairAtomInvGen,
              (python::arg("includeChirality") = false),
              "Get an atom pair atom-invariant generator\n\n"
              "  ARGUMENTS:\n"
              "    - includeChirality: if set, chirality will be taken into "
              "account for invariants\n"
              "  RETURNS: AtomInvariantsGenerator\n\n",
              python::return_value_policy<python::manage_new_object>());

  return;
}
}  // namespace AtomPairWrapper

}  // namespace RDKit
