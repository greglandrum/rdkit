//
//  Copyright (C) 2025 RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_AVALONFPGEN_H_2025
#define RD_AVALONFPGEN_H_2025

#include <GraphMol/Fingerprints/FingerprintGenerator.h>

namespace RDKit {
namespace AvalonFP {

//! Default feature set: all of the Avalon features that do not depend on
//! scaffold information (the equivalent of Avalon's USE_ALL_FEATURES)
const std::uint32_t defaultAvalonBitFlags = 0x7FFF;

/*!
  Computes the Avalon feature counts of a molecule. This is a pure RDKit
  port of the Avalon toolkit's SetFingerprintCountsWithFocus().

  \param mol       the molecule
  \param fpSize    number of slots (hash bins) in the result
  \param bitFlags  controls which features are used (Avalon's \c USE_* flags)
  \param isQuery   hydrogen counts are only taken from explicit hydrogens and
                   the substitution dependent features are disabled
  \param focusAtom if >=0 (an atom index), the counts are computed as if this
                   atom (and its bonds) were not part of the molecule

  \return a vector of length fpSize with the number of features hashed into
  each slot. A bit-vector fingerprint (Avalon's SetFingerprintBits()) is the
  set of slots with a non-zero count.

  \param accumulateAsQuery  for non-query molecules the algorithm is run twice
                   (Avalon's aromaticity model, then Daylight-like
                   aromaticity) and the results are accumulated. This
                   controls whether the second run uses query mode, which is
                   what AvalonTools::getAvalonCountFP() does (but not
                   AvalonTools::getAvalonFP()).

  Aromaticity perception follows the Avalon toolkit, but ring perception and
  implicit hydrogen counts come from RDKit, so results may differ from the
  Avalon toolkit in unusual cases (e.g. cage ring systems, unusual valences).
*/
RDKIT_FINGERPRINTS_EXPORT std::vector<std::uint32_t> getAvalonCounts(
    const ROMol &mol, std::uint32_t fpSize = 512,
    std::uint32_t bitFlags = defaultAvalonBitFlags, bool isQuery = false,
    int focusAtom = -1, bool accumulateAsQuery = false);

class RDKIT_FINGERPRINTS_EXPORT AvalonArguments : public FingerprintArguments {
 public:
  std::uint32_t d_bitFlags = defaultAvalonBitFlags;
  bool df_isQuery = false;
  bool df_accumulateAsQuery = false;

  std::string infoString() const override;
  void toJSON(boost::property_tree::ptree &pt) const override;
  void fromJSON(const boost::property_tree::ptree &pt) override;

  /*!
   \param bitFlags  Avalon \c USE_* feature flags
   \param isQuery   generate the fingerprint for a query molecule
   \param countSimulation  use count simulation (not normally useful here)
   \param countBounds      bounds for count simulation
   \param fpSize    size of the fingerprint
   \param accumulateAsQuery see getAvalonCounts()
  */
  AvalonArguments(std::uint32_t bitFlags = defaultAvalonBitFlags,
                  bool isQuery = false, bool countSimulation = false,
                  const std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
                  std::uint32_t fpSize = 512, bool accumulateAsQuery = false);
};

template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT AvalonAtomEnv
    : public AtomEnvironment<OutputType> {
  const OutputType d_bitId;

 public:
  OutputType getBitId(FingerprintArguments *,
                      const std::vector<std::uint32_t> *,
                      const std::vector<std::uint32_t> *, AdditionalOutput *,
                      bool hashResults = false,
                      const std::uint64_t fpSize = 0) const override;
  void updateAdditionalOutput(AdditionalOutput *,
                              std::uint64_t) const override;
  explicit AvalonAtomEnv(const OutputType bitId) : d_bitId(bitId) {}
};

template <typename OutputType>
class RDKIT_FINGERPRINTS_EXPORT AvalonEnvGenerator
    : public AtomEnvironmentGenerator<OutputType> {
  std::uint32_t d_fpSize = 512;

 public:
  AvalonEnvGenerator() = default;
  explicit AvalonEnvGenerator(std::uint32_t fpSize) : d_fpSize(fpSize) {}

  //! fromAtoms and ignoreAtoms are not supported (an exception is thrown if
  //! they are non-empty); the other arguments are ignored
  std::vector<AtomEnvironment<OutputType> *> getEnvironments(
      const ROMol &mol, FingerprintArguments *arguments,
      const std::vector<std::uint32_t> *fromAtoms,
      const std::vector<std::uint32_t> *ignoreAtoms, int confId,
      const AdditionalOutput *additionalOutput,
      const std::vector<std::uint32_t> *atomInvariants,
      const std::vector<std::uint32_t> *bondInvariants,
      bool hashResults = false) const override;

  std::string infoString() const override;
  void toJSON(boost::property_tree::ptree &pt) const override;
  void fromJSON(const boost::property_tree::ptree &pt) override;

  OutputType getResultSize() const override;
};

/*!
 \brief Get an Avalon fingerprint generator (pure RDKit implementation)

 The bit-vector fingerprints (getFingerprint()) correspond to Avalon's
 SetFingerprintBits(), the count fingerprints (getCountFingerprint() and
 getSparseCountFingerprint()) to SetFingerprintCountsWithFocus().

 \param fpSize    size of the generated fingerprint
 \param bitFlags  Avalon \c USE_* feature flags
 \param isQuery   generate the fingerprint for a query molecule
 \param accumulateAsQuery  see getAvalonCounts(); set this to reproduce the
                  count fingerprints of AvalonTools::getAvalonCountFP()
*/
template <typename OutputType>
RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<OutputType> *getAvalonGenerator(
    std::uint32_t fpSize = 512, std::uint32_t bitFlags = defaultAvalonBitFlags,
    bool isQuery = false, bool accumulateAsQuery = false);
// \overload
template <typename OutputType>
RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<OutputType> *getAvalonGenerator(
    const AvalonArguments &args);

}  // namespace AvalonFP
}  // namespace RDKit

#endif
