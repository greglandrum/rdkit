#include <algorithm>
#include "ExtendedStereoGroup.h"

namespace RDKit {

ExtendedStereoGroup::ExtendedStereoGroup(ExtendedStereoGroupType grouptype,
                                         std::vector<Atom *> &&atoms)
    : d_grouptype(grouptype), d_atoms(atoms) {}
ExtendedStereoGroup::ExtendedStereoGroup(ExtendedStereoGroupType grouptype,
                                         const std::vector<Atom *> &atoms)
    : d_grouptype(grouptype), d_atoms(atoms) {}

ExtendedStereoGroupType ExtendedStereoGroup::getGroupType() const {
  return d_grouptype;
}

const std::vector<Atom *> &ExtendedStereoGroup::getAtoms() const {
  return d_atoms;
}

void removeGroupsWithAtom(const Atom *atom,
                          std::vector<ExtendedStereoGroup> &groups) {
  auto containsAtom = [atom](const ExtendedStereoGroup &group) {
    return std::find(group.getAtoms().cbegin(), group.getAtoms().cend(),
                     atom) != group.getAtoms().cend();
  };
  groups.erase(std::remove_if(groups.begin(), groups.end(), containsAtom),
               groups.end());
}

void removeGroupsWithAtoms(const std::vector<Atom *> &atoms,
                           std::vector<ExtendedStereoGroup> &groups) {
  auto containsAnyAtom = [atoms](const ExtendedStereoGroup &group) {
    for (auto atom : atoms) {
      if (std::find(group.getAtoms().cbegin(), group.getAtoms().cend(), atom) !=
          group.getAtoms().cend()) {
        return true;
      }
    }
    return false;
  };
  groups.erase(std::remove_if(groups.begin(), groups.end(), containsAnyAtom),
               groups.end());
}

}  // namespace RDKit
