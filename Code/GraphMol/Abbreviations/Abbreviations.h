//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_ABBREVIATIONS_H
#define RD_ABBREVIATIONS_H
#include <vector>
#include <string>
#include <memory>

namespace RDKit {
class ROMol;
class RWMol;

namespace Abbreviations {
RDKIT_ABBREVIATIONS_EXPORT struct AbbreviationDefinition {
  std::string llabel;
  std::string rlabel;
  std::string smarts;
  std::shared_ptr<ROMol> mol;
  bool operator==(const AbbreviationDefinition& other) const {
    return llabel == other.llabel && rlabel == other.rlabel &&
           smarts == other.smarts;
  }
};
RDKIT_ABBREVIATIONS_EXPORT struct AbbreviationMatch {
  std::vector<std::pair<int, int>> match;
  AbbreviationDefinition abbrev;
  AbbreviationMatch(const std::vector<std::pair<int, int>>& matchArg,
                    const AbbreviationDefinition& abbrevArg)
      : match(matchArg), abbrev(abbrevArg){};
};
namespace common_properties {
RDKIT_ABBREVIATIONS_EXPORT extern const std::string numDummies;
}
namespace Utils {
//! returns the default set of abbreviation definitions
RDKIT_ABBREVIATIONS_EXPORT std::vector<AbbreviationDefinition>
getDefaultAbbreviations();
//! returns the default set of linker definitions
RDKIT_ABBREVIATIONS_EXPORT std::vector<AbbreviationDefinition>
getDefaultLinkers();

//! parses a string describing abbreviation matches and returns the result
/*

\param text the data to be parsed, see below for the format
\param removeExtraDummies controls whether or not dummy atoms beyond atom 0 are
       removed. Set this to true to create abbreviations for linkers
\param allowConnectionToDummies allows abbreviations to directly connect to
       abbreviations. set this to true for linkers

Format of the text data:
  A series of lines, each of which contains:

    label rlabel SMARTS

  where label is the label used for the abbreviation,
  rlabel is the display label if a bond comes in from the right,
  and SMARTS is the SMARTS definition of the abbreviation.
  Use dummies to indicate attachment points. The assumption is that the first
  atom is a dummy (one will be added if this is not true) and that the second
  atom is the surrogate for the rest of the group.

*/
RDKIT_ABBREVIATIONS_EXPORT std::vector<AbbreviationDefinition>
parseAbbreviations(const std::string& text, bool removeExtraDummies = false,
                   bool allowConnectionToDummies = false);
}  // namespace Utils

//! returns all matches for the abbreviations across the molecule
/*!

    \param abbrevs the abbreviations to look for. This list is used in order.
    \param maxCoverage any abbreviation that covers than more than this fraction
        of the molecule's atoms (not counting dummies) will not be returned.
*/
RDKIT_ABBREVIATIONS_EXPORT std::vector<AbbreviationMatch>
findApplicableAbbreviationMatches(
    const ROMol& mol, const std::vector<AbbreviationDefinition>& abbrevs,
    double maxCoverage = 0.4);
//! applies the abbreviation matches to a molecule, modifying it in place.
RDKIT_ABBREVIATIONS_EXPORT void applyMatches(
    RWMol& mol, const std::vector<AbbreviationMatch>& matches);
//! creates "SUP" SubstanceGroups on the molecule describing the abbreviation
RDKIT_ABBREVIATIONS_EXPORT void labelMatches(
    RWMol& mol, const std::vector<AbbreviationMatch>& matches);
//! convenience function for finding and applying abbreviations
RDKIT_ABBREVIATIONS_EXPORT void condenseMolAbbreviations(
    RWMol& mol, const std::vector<AbbreviationDefinition>& abbrevs,
    double maxCoverage = 0.4, bool sanitize = true);
//! convenience function for finding and labeling abbreviations
RDKIT_ABBREVIATIONS_EXPORT void labelMolAbbreviations(
    RWMol& mol, const std::vector<AbbreviationDefinition>& abbrevs,
    double maxCoverage = 0.4);

}  // namespace Abbreviations
}  // namespace RDKit
#endif
