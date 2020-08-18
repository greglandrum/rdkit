//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY
#include <RDBoost/iterator_next.h>
#include <RDBoost/python.h>

#include <string>

// ours
#include <GraphMol/FileParsers/MultithreadedSDMolSupplier.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/FileParseException.h>

#include "MultithreadedMolSupplier.h"

namespace python = boost::python;

namespace RDKit {

std::string multiSDMolSupplierClassDoc =
    "A class which concurrently supplies molecules from a text file.\n\
\n\
  Usage examples:\n\
\n\
    1) Lazy evaluation: the molecules might not be constructed until we ask for them:\n\n\
       >>> suppl = MultithreadedSDMolSupplier('in.smi')\n\
       >>> while not suppl.atEnd():\n\
			 ...    mol = next(mol)\n\
       ...    mol.GetNumAtoms()\n\
\n\
  If the input file has a title line and more than two columns (smiles and id), the\n\
  additional columns will be used to set properties on each molecule.  The properties\n\
  are accessible using the mol.GetProp(propName) method.\n\
\n";

std::string multiSdsDocStr =
    "Constructor\n\n\
  ARGUMENTS: \n\
\n\
    - fileName: name of the file to be read\n\
\n\
    - sanitize: (optional) toggles sanitization of molecules as they are read.\n\
      Defaults to 1.\n\
\n\
    - removeHs: (optional) removes Hs\n\
\n\
    - strictParsing: (optional) allows strict or lax parsing\n\
\n\
    - numWriterThreads: (optional) number of writer threads\n\
\n\
    - sizeInputQueue: (optional) size of input queue\n\
\n\
    - sizeOutputQueue: (optional) size of output queue\n\
\n";

struct multiSDMolSup_wrap {
  static void wrap() {
    python::class_<MultithreadedSDMolSupplier, boost::noncopyable>(
        "MultithreadedSDMolSupplier", multiSDMolSupplierClassDoc.c_str(),
        python::init<>())
        .def(python::init<std::string, bool, bool, bool, unsigned int, size_t,
                          size_t>(
            (python::arg("fileName"), python::arg("sanitize") = true,
             python::arg("removeHs") = true,
             python::arg("strictParsing") = true,
             python::arg("numWriterThreads") = 2,
             python::arg("sizeInputQueue") = 5,
             python::arg("sizeOutputQueue") = 5),
            multiSdsDocStr.c_str()))
        .def("__next__",
             (ROMol * (*)(MultithreadedSDMolSupplier *)) & MolForwardSupplNext,
             "Returns the next molecule in the file. Raises _StopIteration_ "
             "on EOF.\n",
             python::return_value_policy<python::manage_new_object>())
        .def("atEnd", &MultithreadedSDMolSupplier::atEnd,
             "Returns true if we have read all records else false.\n")
        .def("GetLastRecordId",
             (unsigned int (*)(MultithreadedSDMolSupplier *)) & MolSupplLastId,
             "Returns the record id for the last extracted item.\n")
        .def("GetItemText",
             (std::string(*)(MultithreadedSDMolSupplier *)) & MolSupplLastItem,
             "Returns the text for the last extracted item.\n");
  };
};
}  // namespace RDKit

void wrap_multiSDSupplier() { RDKit::multiSDMolSup_wrap::wrap(); }
