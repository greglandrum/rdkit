//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RDKIT_HYPERSPACESUBSTRUCTUREDETAILS_H
#define RDKIT_HYPERSPACESUBSTRUCTUREDETAILS_H

#include <vector>

namespace RDKit {
class ROMol;
namespace HyperspaceSearch {

class Hyperspace;

namespace details {
// Find all combinations of M things selected from N.
RDKIT_HYPERSPACESEARCH_EXPORT std::vector<std::vector<unsigned int>> combMFromN(
    unsigned int m, unsigned int n);
// Find all permutations of M things selected from N.
RDKIT_HYPERSPACESEARCH_EXPORT std::vector<std::vector<unsigned int>> permMFromN(
    unsigned int m, unsigned int n);
// Split the molecule into pieces breaking up to maxBondSplits bonds at a
// time.
RDKIT_HYPERSPACESEARCH_EXPORT std::vector<std::vector<std::unique_ptr<ROMol>>>
splitMolecule(const ROMol &query, unsigned int maxBondSplits);
// Counts the number of [1*], [2*]...[4*] in the string.
RDKIT_HYPERSPACESEARCH_EXPORT int countConnections(const std::string &smiles);
}  // namespace details
}  // namespace HyperspaceSearch
}  // namespace RDKit

#endif  // RDKIT_HYPERSPACESUBSTRUCTURESEARCH_H
