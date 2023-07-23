//
// Copyright (C) David Cosgrove 2023
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

// A class to hold the results of a RASCAL MCES determination
// between 2 molecules.  Contains the bonds and atoms that
// correspond between the molecules, and also a SMARTS pattern
// defining the MCES.
//

#ifndef RASCALRESULT_H
#define RASCALRESULT_H

#include <vector>

#include <GraphMol/ROMol.h>

namespace RDKit {

namespace RascalMCES {

class RascalResult {
 public:
  RascalResult(const RDKit::ROMol &mol1, const RDKit::ROMol &mol2,
               const std::vector<std::vector<int>> &adjMatrix1,
               const std::vector<std::vector<int>> &adjMatrix2,
               const std::vector<unsigned int> &clique,
               const std::vector<std::pair<int, int>> &vtx_pairs, bool timedOut,
               bool swapped, bool chiralSmarts, bool ringMatchesRingOnly,
               bool singleLargestFrag, int minFragSep);

  RascalResult(const RascalResult &other);

  RascalResult(RascalResult &&other) = default;

  ~RascalResult() = default;

  RascalResult &operator=(const RascalResult &other);

  RascalResult &operator=(RascalResult &&other) = default;

  // Cut the result down to the single largest fragment.  This is
  // irrecoverably destructive.
  void largestFragOnly();

  std::shared_ptr<RDKit::ROMol> mol1() const { return d_mol1; };

  std::shared_ptr<RDKit::ROMol> mol2() const { return d_mol2; };

  std::vector<std::pair<int, int>> bondMatches() const { return d_bondMatches; }

  std::vector<std::pair<int, int>> atomMatches() const { return d_atomMatches; }

  int numFrags() const;

  int ringNonRingBondScore() const;

  int atomMatchScore() const;

  int maxDeltaAtomAtomDist() const;

  // returns number of atoms for largest fragment.
  int largestFragSize() const;

  std::string smarts() const;

  bool timedout() const { return d_timedOut; };

  double similarity() const;

 private:
  std::shared_ptr<RDKit::ROMol> d_mol1;
  std::shared_ptr<RDKit::ROMol> d_mol2;
  std::vector<std::pair<int, int>> d_bondMatches;
  std::vector<std::pair<int, int>> d_atomMatches;

  mutable std::string d_smarts;
  bool d_timedOut{false};
  bool d_chiralSmarts{false};
  bool d_ringMatchesRingOnly{false};
  int d_maxFragSep{-1};

  // These are used for sorting the results.
  mutable int d_numFrags{-1};
  mutable int d_ringNonRingBondScore{-1};
  mutable int d_atomMatchScore{-1};
  mutable int d_maxDeltaAtomAtomDist{-1};
  mutable int d_largestFragSize{-1};

  std::string createSmartsString() const;

  void matchCliqueAtoms(const std::vector<std::vector<int>> &mol1_adj_matrix);

  // If the clique involves a fragment that is more than d_maxFragSep from
  // any other frag in either molecule, discard the smaller frag.
  void applyMaxFragSep();

  // Make the fragments for either mol1 or mol2.  If molNum is not 1 or 2,
  // returns nullptr.
  RDKit::ROMol *makeMolFrags(int molNum) const;

  int calcRingNonRingScore() const;

  int calcAtomMatchScore() const;

  int calcLargestFragSize() const;

  // If there are multiple fragments, can be helpful as a tie-breaker.  It's the
  // maximum difference between through-bond distances between matching atoms in
  // the 2 molecules.
  int calcMaxDeltaAtomAtomDistScore() const;
};

bool resultSort(const RascalResult &res1, const RascalResult &res2);

void extractClique(const std::vector<unsigned int> &clique,
                   const std::vector<std::pair<int, int>> &vtxPairs,
                   bool swapped, std::vector<std::pair<int, int>> &bondMatches);

// do some simple cleaning of the SMARTS, to make it more user-friendly.
void cleanSmarts(std::string &smarts);

// Primarily for debugging, these write out the corresponding bonds/atoms
// in Python list format, for ease of cut/paste into a highlighted image
// creation.
void printBondMatches(const RascalResult &res, std::ostream &os);

void printAtomMatches(const RascalResult &res, std::ostream &os);

// Calculate the Johnson similarity between the two molecules using the given
// bondMatches.  It's the fraction of the 2 molecules that are in common,
// somewhat akin to the tanimoto - the square of the number of atoms plus
// number of bonds in the MCES divided by the product of the sums of the number
// of atoms and bonds in the 2 molecules.
// It has nothing to do with lying UK politicians.
double johnsonSimilarity(const std::vector<std::pair<int, int>> &bondMatches,
                         const RDKit::ROMol &mol1, const RDKit::ROMol &mol2);

}  // namespace RascalMCES
}  // namespace RDKit

#endif  // RASCALRESULT_H
