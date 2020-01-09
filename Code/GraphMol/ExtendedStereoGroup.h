//
//  Copyright (C) 2018 T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file ExtendedStereoGroup.h

  \brief Defines the class ExtendedStereoGroup which stores relationships
  between the absolute configurations of atoms within a structure.

*/

#include <RDGeneral/export.h>
#ifndef RD_ExtendedStereoGroup_092018
#define RD_ExtendedStereoGroup_092018

#include <vector>

namespace RDKit {
class Atom;

// OR means that it is known to be one or the other, but not both
// AND means that it is known to be a mix.
enum class ExtendedStereoGroupType {
  STEREO_ABSOLUTE = 0,
  STEREO_OR = 1,
  STEREO_AND = 2
};

//! ExtendedStereoGroup is a collection of atoms with a known stereochemical
//! relationship
/*!
  Used to help represent a sample with unknown stereochemistry, or that is a mix
  of diastereomers.

 */
class RDKIT_GRAPHMOL_EXPORT ExtendedStereoGroup {
 private:
  ExtendedStereoGroupType d_grouptype;
  std::vector<Atom*> d_atoms;

 public:
  ExtendedStereoGroup()
      : d_grouptype(ExtendedStereoGroupType::STEREO_ABSOLUTE), d_atoms(0u){};
  // Takes control of atoms if possible.
  ExtendedStereoGroup(ExtendedStereoGroupType grouptype,
                      std::vector<Atom*>&& atoms);
  ExtendedStereoGroup(ExtendedStereoGroupType grouptype,
                      const std::vector<Atom*>& atoms);
  ExtendedStereoGroupType getGroupType() const;
  const std::vector<Atom*>& getAtoms() const;
  // Seems odd to have to define these, but otherwise the SWIG wrappers
  // won't build
  bool operator==(const ExtendedStereoGroup& other) const {
    return (d_grouptype == other.d_grouptype) && (d_atoms == other.d_atoms);
  };
  bool operator!=(const ExtendedStereoGroup& other) const {
    return (d_grouptype != other.d_grouptype) || (d_atoms != other.d_atoms);
  };
};
RDKIT_GRAPHMOL_EXPORT void removeGroupsWithAtom(
    const Atom* atom, std::vector<ExtendedStereoGroup>& groups);
RDKIT_GRAPHMOL_EXPORT void removeGroupsWithAtoms(
    const std::vector<Atom*>& atoms, std::vector<ExtendedStereoGroup>& groups);

}  // namespace RDKit

#endif
