#ifndef RD_FINGERPRINTGEN_H_2018_05
#define RD_FINGERPRINTGEN_H_2018_05

#include <DataStructs/SparseIntVect.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseBitVect.h>

namespace RDKit {
class ROMol;

union AdditionalOutput {
  // will hold differently stuctured additonal output

  std::map<boost::uint32_t,
           std::vector<std::pair<boost::uint32_t, boost::uint32_t>>>
      *bitInfoMap;
  // morgan fp
  // maps bitId -> vector of (atomId, radius)

  std::pair<std::vector<std::vector<boost::uint32_t>>,
            std::map<boost::uint32_t, std::vector<std::vector<int>>>> *bitInfo;
  // rdkit fp
  // first part, vector of bits set for each atom, must have the same size as
  // atom count for molecule
  // second part, maps bitId -> vector of paths

  std::vector<unsigned int> *atomCounts;
  // number of paths that set bits for each atom, must have the same size as
  // atom count for molecule
};

/*
 Class that holds the arguments for fingerprint generation
 countSimulation: whether or not to return count simulated fingerprint
*/
class FingerprintArguments {
 protected:
  FingerprintArguments(const bool countSimulation);

 public:
  const bool countSimulation;

  /*
    Returns the size of the fingerprint based on arguments
  */
  virtual unsigned int getResultSize() const = 0;

  virtual ~FingerprintArguments() = 0;
};

/*
Class to hold atom environments and hash them into bit ids
*/
class AtomEnvironment {
 public:
  /*
  returns hashed bit id from this atom environment
  arguments: fingerprint style specific arguments
  atomInvariants: a list of invariants to use for the atom hashes
  bondInvariants: a list of invariants to use for the bond hashes
  */
  virtual boost::uint32_t getBitId(
      FingerprintArguments *arguments,
      const std::vector<boost::uint32_t> *atomInvariants = 0,
      const std::vector<boost::uint32_t> *bondInvariants = 0) const = 0;

  virtual ~AtomEnvironment() = 0;
};

/*
    Class for generating atom environments from the molecule and clean up
   resources after atom environments are no longer needed
*/
class AtomEnvironmentGenerator {
 public:
  /*
  Returns atom environments for a given molecule
  mol: molecule to generate environments from
  arguments: fingerprinting type specific arguments
  fromAtoms: indicies of atoms to be used during fingerprint generation, usage
  depends on fingerprint style
  ignoreAtoms: indicies of atoms to be ignored during fingerprint generation,
  usage depends on fingerprint style
  confId: 3D configuration to use during fingerprint generation
  additionalOutput: used to store additional outputs
  atomInvariants: a list of invariants to use for the atom hashes
  bondInvariants: a list of invariants to use for the bond hashes
  */
  virtual std::vector<AtomEnvironment *> getEnvironments(
      const ROMol &mol, FingerprintArguments *arguments,
      const std::vector<boost::uint32_t> *fromAtoms = 0,
      const std::vector<boost::uint32_t> *ignoreAtoms = 0,
      const int confId = -1, const AdditionalOutput *additionalOutput = 0,
      const std::vector<boost::uint32_t> *atomInvariants = 0,
      const std::vector<boost::uint32_t> *bondInvariants = 0) const = 0;

  /*
  Does clean up for shared resources created for atom environments to use, does
  not remove instance in the vector just the resources they use if any
  atomEnvironments: vector of atom environments that hold the resources to be
  cleaned
  */
  virtual void cleanUpEnvironments(
      std::vector<AtomEnvironment *> atomEnvironments) const = 0;

  virtual ~AtomEnvironmentGenerator() = 0;
};

class AtomInvariantsGenerator {
  // arguments

 public:
  virtual std::vector<boost::uint32_t> getAtomInvariants(
      const ROMol &mol) const = 0;
};

class BondInvariantsGenerator {
  // arguments

 public:
  virtual std::vector<boost::uint32_t> getBondInvariants(
      const ROMol &mol) const = 0;
};

/*
Class to generate same fingerprint style for different output formats
atomEnvironmentGenerator: environment generator class to use for fingerprinting
atomInvariantsGenerator: atom invariants generator
bondInvariantsGenerator: bond invariants generator
fingerprintArguments: holds fingerprint style specific arguments
*/
class FingerprintGenerator {
  FingerprintArguments *fingerprintArguments;
  AtomEnvironmentGenerator *atomEnvironmentGenerator;
  AtomInvariantsGenerator *atomInvariantsGenerator;
  BondInvariantsGenerator *bondInvariantsGenerator;

 public:
  FingerprintGenerator(AtomEnvironmentGenerator *atomEnvironmentGenerator,
                       FingerprintArguments *fingerprintArguments,
                       AtomInvariantsGenerator *atomInvariantsGenerator = 0,
                       BondInvariantsGenerator *bondInvariantsGenerator = 0);

  /*
  Deletes atomEnvironmentGenerator, fingerprintArguments
  since they have virtual destructors they will be properly deleted but the
  instance can not be used anymore after calling this
  */
  void cleanUpResources();

  SparseIntVect<boost::uint32_t> *getFingerprint(
      const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms = 0,
      const std::vector<boost::uint32_t> *ignoreAtoms = 0,
      const int confId = -1,
      const AdditionalOutput *additionalOutput = 0) const;

  SparseBitVect *getFingerprintAsBitVect(
      const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms = 0,
      const std::vector<boost::uint32_t> *ignoreAtoms = 0,
      const int confId = -1,
      const AdditionalOutput *additionalOutput = 0) const;

  SparseIntVect<boost::uint32_t> *getFoldedFingerprint(
      const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms = 0,
      const std::vector<boost::uint32_t> *ignoreAtoms = 0,
      const AdditionalOutput *additionalOutput = 0) const;

  ExplicitBitVect *getFoldedFingerprintAsBitVect(
      const ROMol &mol, const std::vector<boost::uint32_t> *fromAtoms = 0,
      const std::vector<boost::uint32_t> *ignoreAtoms = 0,
      const AdditionalOutput *additionalOutput = 0) const;
};
}  // namespace RDKit

#endif