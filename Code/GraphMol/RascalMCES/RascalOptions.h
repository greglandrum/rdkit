//
// Copyright (C) David Cosgrove 2023
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RASCALOPTIONS_H
#define RASCALOPTIONS_H

namespace RDKit {

namespace RascalMCES {

struct RascalOptions {
  double similarityThreshold =
      0.7;    // if calculated below this, no MCES will be evaluated.
  bool completeAromaticRings =
      true;   // if true, partial aromatic rings won't be returned
  bool ringMatchesRingOnly =
      false;  // if true, ring bonds won't match non-ring bonds
  bool exactChirality = false;  // if true, R must match R and S match S.
  bool singleLargestFrag =
      false; /* if true, only return a single fragment for the MCES. Default
                is to produce multiple matching fragments if necessary. */
  int minFragSize =
      -1;    /* minimum number of atoms in any fragment - -1 means no minimum */
  int maxFragSeparation = -1; /* biggest through-bond distance that bonds can
                               match. -1 means no limit. */
  bool allBestMCESs =
      false; /* If true, all MCESs are returned, in order of diminishing score.
                This is likely to result in higher run times. */
  int timeout = 60;  // max run time, in seconds. -1 means no max.
  bool doEquivBondPruning =
      false;         /* This might make the code run a bit faster in some
                        circumstances, but on average it is very marginal. */
};
}  // namespace RascalMCES
}  // namespace RDKit

#endif  // RASCALOPTIONS_H
