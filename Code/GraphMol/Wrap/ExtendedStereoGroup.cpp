//
//  Copyright (C) 2018 Dan Nealschneider
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY
#include <boost/python.hpp>
#include <string>

// ours
#include <RDBoost/Wrap.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ExtendedStereoGroup.h>

namespace python = boost::python;

namespace RDKit {

namespace {
std::string ExtendedStereoGroupClassDoc =
    "A collection of atoms with a defined stereochemical relationship.\n\n"
    "Used to help represent a sample with unknown stereochemistry, or that "
    "is a mix\nof diastereomers.\n";

ExtendedStereoGroup *createExtendedStereoGroup(ExtendedStereoGroupType typ, ROMol &mol,
                               python::object atomIds) {
  std::vector<Atom *> cppAtoms;
  python::stl_input_iterator<unsigned int> beg(atomIds), end;
  while (beg != end) {
    unsigned int v = *beg;
    if (v >= mol.getNumAtoms())
      throw_value_error("atom index exceeds mol.GetNumAtoms()");
    cppAtoms.push_back(mol.getAtomWithIdx(v));
    ++beg;
  }
  ExtendedStereoGroup *sg = new ExtendedStereoGroup(typ, cppAtoms);
  return sg;
}

python::object getAtomsHelper(ExtendedStereoGroup &sg) {
  python::list res;
  for (auto at : sg.getAtoms()) {
    res.append(boost::ref(*at));
  }
  return python::tuple(res);
}
}  // namespace

struct ExtendedStereoGroup_wrap {
  static void wrap() {
    python::enum_<RDKit::ExtendedStereoGroupType>("ExtendedStereoGroupType")
        .value("STEREO_ABSOLUTE", RDKit::ExtendedStereoGroupType::STEREO_ABSOLUTE)
        .value("STEREO_OR", RDKit::ExtendedStereoGroupType::STEREO_OR)
        .value("STEREO_AND", RDKit::ExtendedStereoGroupType::STEREO_AND)
        .export_values();

    python::class_<ExtendedStereoGroup, boost::shared_ptr<ExtendedStereoGroup>>(
        "ExtendedStereoGroup", ExtendedStereoGroupClassDoc.c_str(), python::no_init)
        .def("GetGroupType", &ExtendedStereoGroup::getGroupType,
             "Returns the ExtendedStereoGroupType.\n")
        .def("GetAtoms", getAtomsHelper,
             "access the atoms in the ExtendedStereoGroup.\n");

    python::def("CreateExtendedStereoGroup", &createExtendedStereoGroup,
                "creates a ExtendedStereoGroup associated with a molecule from a list "
                "of atom Ids",
                (python::arg("ExtendedStereoGroupType"), python::arg("mol"),
                 python::arg("atomIds")),
                python::return_value_policy<
                    python::manage_new_object,
                    python::with_custodian_and_ward_postcall<0, 2>>());
  }
};
}  // namespace RDKit

void wrap_ExtendedStereoGroup() { RDKit::ExtendedStereoGroup_wrap::wrap(); }
