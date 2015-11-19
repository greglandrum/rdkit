//
//  Copyright (C) 2015 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include "FileParsers.h"
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>

#include <string>

using namespace RDKit;

namespace local_data {
// basic example from Brian Cole
std::string basic_json_ethane =
    "{"
    "   \"type\":\"mol\","
    "   \"version\":\"0.9\","
    "   \"title\":\"ethane\","
    "   \"dimension\":2,"
    "   \"atomdefaults\":{"
    "      \"implicithcount\":3,"
    "      \"formalcharge\":0,"
    "      \"stereo\":\"Undefined\""
    "   },"
    "   \"atoms\":["
    "      {"
    "         \"element\":6,"
    "         \"coords\":["
    "            0.0,"
    "            0.0"
    "         ]"
    "      },"
    "      {"
    "         \"element\":6,"
    "         \"coords\":["
    "            0.0,"
    "            1.0"
    "         ]"
    "      }"
    "   ],"
    "   \"bonddefaults\":{"
    "      \"stereo\" : \"Undefined\","
    "      \"order\":1"
    "   },"
    "   \"bonds\":["
    "      {"
    "         \"atoms\":["
    "            0,"
    "            1"
    "         ],"
    "         \"order\":1"
    "      }"
    "   ]"
    "   }";
// basic example from Brian Cole
std::string basic_json_phenol =
    "{"
    "   \"type\":\"mol\","
    "   \"version\":\"0.9\","
    "   \"title\":\"phenol\","
    "   \"dimension\":2,"
    "   \"atomdefaults\":{"
    "      \"implicithcount\":1,"
    "      \"formalcharge\":0,"
    "      \"stereo\":\"Undefined\""
    "   },"
    "   \"atoms\":["
    "      {"
    "         \"element\":6,"
    "         \"coords\":["
    "            1.735,"
    "            0.0"
    "         ]"
    "      },"
    "      {"
    "         \"element\":6,"
    "         \"coords\":["
    "            0.868,"
    "            -0.498"
    "         ]"
    "      },"
    "      {"
    "         \"element\":6,"
    "         \"coords\":["
    "            0.0,"
    "            0.0"
    "         ]"
    "      },"
    "      {"
    "         \"element\":6,"
    "         \"coords\":["
    "            0.0,"
    "            1.005"
    "         ]"
    "      },"
    "      {"
    "         \"element\":6,"
    "         \"coords\":["
    "            0.867,"
    "            1.513"
    "         ]"
    "      },"
    "      {"
    "         \"element\":6,"
    "         \"implicithcount\":0,"
    "         \"coords\":["
    "            1.735,"
    "            1.005"
    "         ]"
    "      },"
    "      {"
    "         \"element\":8,"
    "         \"coords\":["
    "            2.602,"
    "            1.503"
    "         ]"
    "      }"
    "   ],"
    "   \"bonddefaults\":{"
    "      \"stereo\" : \"Undefined\","
    "      \"order\":1"
    "   },"
    "   \"bonds\":["
    "      {"
    "         \"atoms\":["
    "            0,"
    "            5"
    "         ],"
    "         \"order\":2"
    "      },"
    "      {"
    "         \"atoms\":["
    "            0,"
    "            1"
    "         ]"
    "      },"
    "      {"
    "         \"atoms\":["
    "            1,"
    "            2"
    "         ],"
    "         \"order\":2"
    "      },"
    "      {"
    "         \"atoms\":["
    "            2,"
    "            3"
    "         ]"
    "      },"
    "      {"
    "         \"atoms\":["
    "            3,"
    "            4"
    "         ],"
    "         \"order\":2"
    "      },"
    "      {"
    "         \"atoms\":["
    "            4,"
    "            5"
    "         ]"
    "      },"
    "      {"
    "         \"atoms\":["
    "            5,"
    "            6"
    "         ]"
    "      }"
    "   ]"
    "   }";
}

void testBasics() {
  BOOST_LOG(rdInfoLog) << "testing JSON basics" << std::endl;

  {
    // std::cerr << local_data::basic_json_phenol.substr(1520, 50) << std::endl;
    RWMol *mol = JSONToMol(local_data::basic_json_ethane);
    TEST_ASSERT(mol);
    std::cerr << MolToSmiles(*mol) << std::endl;
    TEST_ASSERT(mol->hasProp("_Name"));
    TEST_ASSERT(mol->getProp<std::string>("_Name") == "ethane");
    TEST_ASSERT(mol->getNumAtoms() == 2);
    TEST_ASSERT(mol->getNumBonds() == 1);
    TEST_ASSERT(mol->getNumConformers() == 1);
    TEST_ASSERT(feq(mol->getConformer().getAtomPos(0).x, 0.0));
    TEST_ASSERT(feq(mol->getConformer().getAtomPos(0).y, 0.0));
    TEST_ASSERT(feq(mol->getConformer().getAtomPos(1).x, 0.0));
    TEST_ASSERT(feq(mol->getConformer().getAtomPos(1).y, 1.0));
    TEST_ASSERT(mol->getAtomWithIdx(0)->getNumExplicitHs() == 3);
    delete mol;
  }

  {
    // std::cerr << local_data::basic_json_phenol.substr(1520, 50) << std::endl;
    RWMol *mol = JSONToMol(local_data::basic_json_phenol);
    TEST_ASSERT(mol);
    std::cerr << MolToSmiles(*mol) << std::endl;
    TEST_ASSERT(mol->hasProp("_Name"));
    TEST_ASSERT(mol->getProp<std::string>("_Name") == "phenol");
    TEST_ASSERT(mol->getNumAtoms() == 7);
    TEST_ASSERT(mol->getNumBonds() == 7);
    TEST_ASSERT(mol->getAtomWithIdx(5)->getNumExplicitHs() == 0);
    TEST_ASSERT(mol->getAtomWithIdx(6)->getNumExplicitHs() == 1);
    delete mol;
  }

  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();

  std::string rdbase = getenv("RDBASE");

  testBasics();

  return 0;
}
