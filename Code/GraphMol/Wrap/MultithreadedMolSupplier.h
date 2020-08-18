//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_WRAP_MOLSUPPLIER_H
#define RD_WRAP_MOLSUPPLIER_H
//! Template functions for wrapping suppliers as python iterators.

#include <GraphMol/RDKitBase.h>
#include <RDBoost/python.h>
#include <RDGeneral/FileParseException.h>

namespace RDKit {
// Note that this returns a pointer to the supplier itself, so be careful
// that it doesn't get deleted by python!
template <typename T>
ROMol *MolForwardSupplNext(T *suppl) {
  ROMol *res = nullptr;
  if (!suppl->atEnd()) {
    res = suppl->next();
  } else {
    PyErr_SetString(PyExc_StopIteration, "End of supplier hit");
    throw boost::python::error_already_set();
  }
  return res;
}

template <typename T>
std::string MolSupplLastItem(T *supp) {
  return supp->getLastItemText();
}

template <typename T>
unsigned int MolSupplLastId(T *supp) {
  return supp->getLastRecordId();
}

}  // namespace RDKit
#endif
