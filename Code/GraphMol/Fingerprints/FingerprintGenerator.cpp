//
//  Copyright (C) 2018 Boran Adas
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <DataStructs/SparseIntVect.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/SparseBitVect.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <RDGeneral/hash/hash.hpp>
#include <cstdint>

namespace RDKit {

template <typename OutputType>
FingerprintArguments<OutputType>::FingerprintArguments(
    const bool countSimulation, const std::vector<std::uint32_t> countBounds,
    std::uint32_t fpSize)
    : d_countSimulation(countSimulation),
      d_countBounds(countBounds),
      d_fpSize(fpSize) {
  PRECONDITION(!countSimulation || !countBounds.empty(),
               "bad count bounds provided");
}

template FingerprintArguments<std::uint32_t>::FingerprintArguments(
    const bool countSimulation, const std::vector<std::uint32_t> countBounds,
    std::uint32_t fpSize);

template FingerprintArguments<std::uint64_t>::FingerprintArguments(
    const bool countSimulation, const std::vector<std::uint32_t> countBounds,
    std::uint32_t fpSize);

template <typename OutputType>
std::string FingerprintArguments<OutputType>::commonArgumentsString() const {
  return "Common arguments\tcountSimulation=" +
         std::to_string(d_countSimulation) +
         " fpSize=" + std::to_string(d_fpSize);
}

template <typename OutputType>
FingerprintGenerator<OutputType>::FingerprintGenerator(
    AtomEnvironmentGenerator<OutputType> *atomEnvironmentGenerator,
    FingerprintArguments<OutputType> *fingerprintArguments,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator, bool ownsAtomInvGenerator,
    bool ownsBondInvGenerator)
    : df_ownsAtomInvGenerator(ownsAtomInvGenerator),
      df_ownsBondInvGenerator(ownsBondInvGenerator) {
  this->dp_atomEnvironmentGenerator = atomEnvironmentGenerator;
  this->dp_fingerprintArguments = fingerprintArguments;
  this->dp_atomInvariantsGenerator = atomInvariantsGenerator;
  this->dp_bondInvariantsGenerator = bondInvariantsGenerator;
}

template FingerprintGenerator<std::uint32_t>::FingerprintGenerator(
    AtomEnvironmentGenerator<std::uint32_t> *atomEnvironmentGenerator,
    FingerprintArguments<std::uint32_t> *fingerprintArguments,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator, bool ownsAtomInvGenerator,
    bool ownsBondInvGenerator);

template FingerprintGenerator<std::uint64_t>::FingerprintGenerator(
    AtomEnvironmentGenerator<std::uint64_t> *atomEnvironmentGenerator,
    FingerprintArguments<std::uint64_t> *fingerprintArguments,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator, bool ownsAtomInvGenerator,
    bool ownsBondInvGenerator);

template <typename OutputType>
FingerprintGenerator<OutputType>::~FingerprintGenerator() {
  delete dp_atomEnvironmentGenerator;
  delete dp_fingerprintArguments;
  if (df_ownsAtomInvGenerator) {
    delete dp_atomInvariantsGenerator;
  }
  if (df_ownsBondInvGenerator) {
    delete dp_bondInvariantsGenerator;
  }
}

template FingerprintGenerator<std::uint32_t>::~FingerprintGenerator();

template FingerprintGenerator<std::uint64_t>::~FingerprintGenerator();

template <typename OutputType>
std::string FingerprintGenerator<OutputType>::infoString() const {
  std::string seperator = " : ";
  return dp_fingerprintArguments->commonArgumentsString() + seperator +
         dp_fingerprintArguments->infoString() + seperator +
         dp_atomEnvironmentGenerator->infoString() + seperator +
         (dp_atomInvariantsGenerator
              ? (dp_atomInvariantsGenerator->infoString() + seperator)
              : "No atom invariants generator") +
         (dp_bondInvariantsGenerator
              ? (dp_bondInvariantsGenerator->infoString() + seperator)
              : "No bond invariants generator");
}

template <typename OutputType>
SparseIntVect<OutputType>
    *FingerprintGenerator<OutputType>::getSparseCountFingerprint(
        const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
        const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
        const AdditionalOutput *additionalOutput,
        const std::vector<std::uint32_t> *customAtomInvariants,
        const std::vector<std::uint32_t> *customBondInvariants) const {
  // create the atom and bond invariants using the generators if there are any,
  // created invariants will be passed to each atom environment's getBitId call
  std::vector<std::uint32_t> *atomInvariants = nullptr;
  if (customAtomInvariants) {
    atomInvariants = new std::vector<std::uint32_t>(*customAtomInvariants);
  } else if (dp_atomInvariantsGenerator) {
    atomInvariants = dp_atomInvariantsGenerator->getAtomInvariants(mol);
  }

  std::vector<std::uint32_t> *bondInvariants = nullptr;
  if (customBondInvariants) {
    bondInvariants = new std::vector<std::uint32_t>(*customBondInvariants);
  } else if (dp_bondInvariantsGenerator) {
    bondInvariants = dp_bondInvariantsGenerator->getBondInvariants(mol);
  }

  // create all atom environments that will generate the bit-ids that will make
  // up the fingerprint
  std::vector<AtomEnvironment<OutputType> *> atomEnvironments =
      dp_atomEnvironmentGenerator->getEnvironments(
          mol, dp_fingerprintArguments, fromAtoms, ignoreAtoms, confId,
          additionalOutput, atomInvariants, bondInvariants);

  // allocate the result
  SparseIntVect<OutputType> *res =
      new SparseIntVect<OutputType>(dp_fingerprintArguments->getResultSize());

  // iterate over every atom environment and generate bit-ids that will make up
  // the fingerprint
  for (auto it = atomEnvironments.begin(); it != atomEnvironments.end(); it++) {
    OutputType bitId = (*it)->getBitId(dp_fingerprintArguments, atomInvariants,
                                       bondInvariants, additionalOutput);

    res->setVal(bitId, res->getVal(bitId) + 1);

    delete (*it);
  }

  delete atomInvariants;
  delete bondInvariants;

  return res;
}

template <typename OutputType>
SparseBitVect *FingerprintGenerator<OutputType>::getSparseFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *customAtomInvariants,
    const std::vector<std::uint32_t> *customBondInvariants) const {
  // generate fingerprint using getSparseCountFingerprint and convert it to
  // SparseBitVect
  SparseIntVect<OutputType> *tempResult = getSparseCountFingerprint(
      mol, fromAtoms, ignoreAtoms, confId, additionalOutput,
      customAtomInvariants, customBondInvariants);
  OutputType countBitsPerBit = dp_fingerprintArguments->d_countBounds.size();

  SparseBitVect *result;

  if (dp_fingerprintArguments->d_countSimulation) {
    OutputType sizeWithCount =
        dp_fingerprintArguments->getResultSize() * countBitsPerBit;
    result = new SparseBitVect(sizeWithCount);
  } else {
    result = new SparseBitVect(dp_fingerprintArguments->getResultSize());
  }

  BOOST_FOREACH (auto val, tempResult->getNonzeroElements()) {
    if (dp_fingerprintArguments->d_countSimulation) {
      for (unsigned int i = 0; i < countBitsPerBit; ++i) {
        // for every bound in the d_countBounds in dp_fingerprintArguments, set
        // a bit if the occurrence count is equal or higher than the bound for
        // that bit
        if (val.second >= dp_fingerprintArguments->d_countBounds[i]) {
          result->setBit(val.first * countBitsPerBit + i);
        }
      }
    } else {
      result->setBit(val.first);
    }
  }

  delete tempResult;
  return result;
}

template <typename OutputType>
SparseIntVect<std::uint32_t>
    *FingerprintGenerator<OutputType>::getCountFingerprint(
        const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
        const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
        const AdditionalOutput *additionalOutput,
        const std::vector<std::uint32_t> *customAtomInvariants,
        const std::vector<std::uint32_t> *customBondInvariants) const {
  // create the atom and bond invariants using the generators if there are any,
  // created invariants will be passed to each atom environment's getBitId call
  std::vector<std::uint32_t> *atomInvariants = nullptr;
  if (customAtomInvariants) {
    atomInvariants = new std::vector<std::uint32_t>(*customAtomInvariants);
  } else if (dp_atomInvariantsGenerator) {
    atomInvariants = dp_atomInvariantsGenerator->getAtomInvariants(mol);
  }

  std::vector<std::uint32_t> *bondInvariants = nullptr;
  if (customBondInvariants) {
    bondInvariants = new std::vector<std::uint32_t>(*customBondInvariants);
  } else if (dp_bondInvariantsGenerator) {
    bondInvariants = dp_bondInvariantsGenerator->getBondInvariants(mol);
  }

  // create all atom environments that will generate the bit-ids that will make
  // up the fingerprint
  std::vector<AtomEnvironment<OutputType> *> atomEnvironments =
      dp_atomEnvironmentGenerator->getEnvironments(
          mol, dp_fingerprintArguments, fromAtoms, ignoreAtoms, confId,
          additionalOutput, atomInvariants, bondInvariants, true);

  // allocate the result
  SparseIntVect<std::uint32_t> *res =
      new SparseIntVect<std::uint32_t>(dp_fingerprintArguments->d_fpSize);

  // iterate over every atom environment and generate bit-ids that will make up
  // the fingerprint
  for (auto it = atomEnvironments.begin(); it != atomEnvironments.end(); it++) {
    OutputType bitId = (*it)->getBitId(dp_fingerprintArguments, atomInvariants,
                                       bondInvariants, additionalOutput, true);

    // dp_fingerprintArguments->d_fpSize is 32 bits so the result will fit
    // into 32 bit integer
    std::uint32_t bitId32 = bitId % dp_fingerprintArguments->d_fpSize;

    res->setVal(bitId32, res->getVal(bitId32) + 1);

    delete (*it);
  }

  delete atomInvariants;
  delete bondInvariants;

  return res;
}

template <typename OutputType>
ExplicitBitVect *FingerprintGenerator<OutputType>::getFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *customAtomInvariants,
    const std::vector<std::uint32_t> *customBondInvariants) const {
  // same idea from getSparseFingerprint, generate the fingerprint using
  // getCountFingerprint and convert it to a ExplicitBitVect
  SparseIntVect<std::uint32_t> *tempResult =
      getCountFingerprint(mol, fromAtoms, ignoreAtoms, confId, additionalOutput,
                          customAtomInvariants, customBondInvariants);
  OutputType countBitsPerBit = dp_fingerprintArguments->d_countBounds.size();

  ExplicitBitVect *result;

  if (dp_fingerprintArguments->d_countSimulation) {
    OutputType sizeWithCount =
        dp_fingerprintArguments->d_fpSize * countBitsPerBit;
    result = new ExplicitBitVect(sizeWithCount);
  } else {
    result = new ExplicitBitVect(dp_fingerprintArguments->d_fpSize);
  }

  BOOST_FOREACH (auto val, tempResult->getNonzeroElements()) {
    // same count simulation logic used in getSparseFingerprint
    if (dp_fingerprintArguments->d_countSimulation) {
      for (unsigned int i = 0; i < countBitsPerBit; ++i) {
        if (val.second >= dp_fingerprintArguments->d_countBounds[i]) {
          result->setBit(val.first * countBitsPerBit + i);
        }
      }
    } else {
      result->setBit(val.first);
    }
  }

  delete tempResult;
  return result;
}

template SparseIntVect<std::uint32_t>
    *FingerprintGenerator<std::uint32_t>::getSparseCountFingerprint(
        const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
        const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
        const AdditionalOutput *additionalOutput,
        const std::vector<std::uint32_t> *customAtomInvariants,
        const std::vector<std::uint32_t> *customBondInvariants) const;

template SparseIntVect<std::uint64_t>
    *FingerprintGenerator<std::uint64_t>::getSparseCountFingerprint(
        const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
        const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
        const AdditionalOutput *additionalOutput,
        const std::vector<std::uint32_t> *customAtomInvariants,
        const std::vector<std::uint32_t> *customBondInvariants) const;

template SparseBitVect *
FingerprintGenerator<std::uint32_t>::getSparseFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *customAtomInvariants,
    const std::vector<std::uint32_t> *customBondInvariants) const;

template SparseBitVect *
FingerprintGenerator<std::uint64_t>::getSparseFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *customAtomInvariants,
    const std::vector<std::uint32_t> *customBondInvariants) const;

template SparseIntVect<std::uint32_t>
    *FingerprintGenerator<std::uint32_t>::getCountFingerprint(
        const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
        const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
        const AdditionalOutput *additionalOutput,
        const std::vector<std::uint32_t> *customAtomInvariants,
        const std::vector<std::uint32_t> *customBondInvariants) const;

template SparseIntVect<std::uint32_t>
    *FingerprintGenerator<std::uint64_t>::getCountFingerprint(
        const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
        const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
        const AdditionalOutput *additionalOutput,
        const std::vector<std::uint32_t> *customAtomInvariants,
        const std::vector<std::uint32_t> *customBondInvariants) const;

template ExplicitBitVect *FingerprintGenerator<std::uint32_t>::getFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *customAtomInvariants,
    const std::vector<std::uint32_t> *customBondInvariants) const;

template ExplicitBitVect *FingerprintGenerator<std::uint64_t>::getFingerprint(
    const ROMol &mol, const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    const AdditionalOutput *additionalOutput,
    const std::vector<std::uint32_t> *customAtomInvariants,
    const std::vector<std::uint32_t> *customBondInvariants) const;

}  // namespace RDKit
