#include <algorithm>
#include "StereoGroup.h"

namespace RDKit {

StereoGroup::StereoGroup(StereoValType grouptype, std::vector<Atom *> &&atoms)
    : d_grouptype(grouptype), d_atoms(atoms) {}
StereoGroup::StereoGroup(StereoValType grouptype,
                         const std::vector<Atom *> &atoms)
    : d_grouptype(grouptype), d_atoms(atoms) {}

StereoValType StereoGroup::getStereoVal() const { return d_grouptype; }

const std::vector<Atom *> &StereoGroup::getAtoms() const { return d_atoms; }

void removeGroupsWithAtom(const Atom *atom, std::vector<StereoGroup> &groups) {
  auto containsAtom = [atom](const StereoGroup &group) {
    return std::find(group.getAtoms().cbegin(), group.getAtoms().cend(),
                     atom) != group.getAtoms().cend();
  };
  groups.erase(std::remove_if(groups.begin(), groups.end(), containsAtom),
               groups.end());
}

void removeGroupsWithAtoms(const std::vector<Atom *> &atoms,
                           std::vector<StereoGroup> &groups) {
  auto containsAnyAtom = [atoms](const StereoGroup &group) {
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
