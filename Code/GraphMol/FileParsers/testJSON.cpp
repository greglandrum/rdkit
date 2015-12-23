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
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>

#include <string>
#include <iostream>
#include <fstream>

using namespace RDKit;

namespace local_data {
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
    "         ]"
    "      }"
    "   ]"
    "   }";
std::string basic_json_ethene =
    "{"
    "   \"type\":\"mol\","
    "   \"version\":\"0.9\","
    "   \"title\":\"ethene\","
    "   \"dimension\":2,"
    "   \"atomdefaults\":{"
    "      \"implicithcount\":2,"
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
    "      \"order\":2"
    "   },"
    "   \"bonds\":["
    "      {"
    "         \"atoms\":["
    "            0,"
    "            1"
    "         ]"
    "      }"
    "   ]"
    "   }";
std::string basic_json_ethanol_anion =
    "{"
    "   \"type\":\"mol\","
    "   \"version\":\"0.9\","
    "   \"title\":\"ethanol\","
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
    "         \"element\":8,"
    "         \"coords\":["
    "            0.0,"
    "            1.0"
    "         ],"
    "         \"implicithcount\":0,"
    "         \"formalcharge\":-1"
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
    "         ]"
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

void testReadBasics() {
  BOOST_LOG(rdInfoLog) << "testing JSON read basics" << std::endl;

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
    RWMol *mol = JSONToMol(local_data::basic_json_ethanol_anion);
    TEST_ASSERT(mol);
    std::cerr << MolToSmiles(*mol) << std::endl;
    TEST_ASSERT(mol->hasProp("_Name"));
    TEST_ASSERT(mol->getProp<std::string>("_Name") == "ethanol");
    TEST_ASSERT(mol->getNumAtoms() == 2);
    TEST_ASSERT(mol->getNumBonds() == 1);
    TEST_ASSERT(mol->getNumConformers() == 1);
    TEST_ASSERT(feq(mol->getConformer().getAtomPos(0).x, 0.0));
    TEST_ASSERT(feq(mol->getConformer().getAtomPos(0).y, 0.0));
    TEST_ASSERT(feq(mol->getConformer().getAtomPos(1).x, 0.0));
    TEST_ASSERT(feq(mol->getConformer().getAtomPos(1).y, 1.0));
    TEST_ASSERT(mol->getAtomWithIdx(0)->getNumExplicitHs() == 3);
    TEST_ASSERT(mol->getAtomWithIdx(1)->getNumExplicitHs() == 0);
    TEST_ASSERT(mol->getAtomWithIdx(0)->getFormalCharge() == 0);
    TEST_ASSERT(mol->getAtomWithIdx(1)->getFormalCharge() == -1);
    TEST_ASSERT(mol->getBondBetweenAtoms(0, 1));
    TEST_ASSERT(mol->getBondBetweenAtoms(0, 1)->getBondType() == Bond::SINGLE);
    delete mol;
  }

  {
    // std::cerr << local_data::basic_json_phenol.substr(1520, 50) << std::endl;
    RWMol *mol = JSONToMol(local_data::basic_json_ethene);
    TEST_ASSERT(mol);
    std::cerr << MolToSmiles(*mol) << std::endl;
    TEST_ASSERT(mol->hasProp("_Name"));
    TEST_ASSERT(mol->getProp<std::string>("_Name") == "ethene");
    TEST_ASSERT(mol->getNumAtoms() == 2);
    TEST_ASSERT(mol->getNumBonds() == 1);
    TEST_ASSERT(mol->getNumConformers() == 1);
    TEST_ASSERT(feq(mol->getConformer().getAtomPos(0).x, 0.0));
    TEST_ASSERT(feq(mol->getConformer().getAtomPos(0).y, 0.0));
    TEST_ASSERT(feq(mol->getConformer().getAtomPos(1).x, 0.0));
    TEST_ASSERT(feq(mol->getConformer().getAtomPos(1).y, 1.0));
    TEST_ASSERT(mol->getAtomWithIdx(0)->getNumExplicitHs() == 2);
    TEST_ASSERT(mol->getBondBetweenAtoms(0, 1));
    TEST_ASSERT(mol->getBondBetweenAtoms(0, 1)->getBondType() == Bond::DOUBLE);
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

void testWriteBasics() {
  BOOST_LOG(rdInfoLog) << "testing JSON write basics" << std::endl;

  {
    // std::cerr << local_data::basic_json_phenol.substr(1520, 50) << std::endl;
    RWMol *mol = SmilesToMol("C1CC1O");
    TEST_ASSERT(mol);
    mol->setProp("_Name", "ethanol");

    std::string json = MolToJSON(*mol);
    std::cerr << json << std::endl;
    delete mol;
  }
  {
    // std::cerr << local_data::basic_json_phenol.substr(1520, 50) << std::endl;
    RWMol *mol = SmilesToMol("c1ccccc1O");
    TEST_ASSERT(mol);
    mol->setProp("_Name", "phenol");
    std::string json = MolToJSON(*mol);
    std::cerr << json << std::endl;
    delete mol;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

namespace {
void _validateFromSmiles(std::string &smiles) {
  RWMol *mol = SmilesToMol(smiles);
  TEST_ASSERT(mol);
  mol->setProp(common_properties::_Name, "validate");
  std::string json = MolToJSON(*mol);
  RWMol *mol2 = JSONToMol(json);
  TEST_ASSERT(mol2);
  TEST_ASSERT(mol2->getProp<std::string>(common_properties::_Name) ==
              "validate");
  TEST_ASSERT(mol2->getNumAtoms() == mol->getNumAtoms());
  TEST_ASSERT(mol2->getNumBonds() == mol->getNumBonds());
  std::string smi1 = MolToSmiles(*mol, true);
  std::string smi2 = MolToSmiles(*mol2, true);
  if (smi1 != smi2) {
    std::cerr << "   " << smi1 << std::endl;
    std::cerr << "   " << smi2 << std::endl;
    mol->debugMol(std::cerr);
    mol2->debugMol(std::cerr);
  }
  TEST_ASSERT(smi1 == smi2);
  delete mol;
  delete mol2;
}
}
void testRoundTripSmiles() {
  BOOST_LOG(rdInfoLog) << "testing JSON round-tripping from Smiles"
                       << std::endl;

  {
    std::string smis[] = {"CCCC",     "C1=CCC1",    "CC[O-]", "C[CH2]",
                          "c1ccccc1", "c1ccc[nH]1", "EOS"};
    for (unsigned int i = 0; smis[i] != "EOS"; ++i) {
      _validateFromSmiles(smis[i]);
    }
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

void testChiralBasics() {
  BOOST_LOG(rdInfoLog) << "testing basic chirality" << std::endl;

  {
    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/basic_chiral1a.json";
    std::ifstream ifs(fName.c_str());
    // std::string json;
    // std::getline(ifs, json);
    // std::cerr << " !!!! " << json << std::endl;

    RWMol *mol = JSONDataStreamToMol(ifs);
    TEST_ASSERT(mol);
    TEST_ASSERT(!mol->hasProp("_Name"));

    std::string smi = MolToSmiles(*mol, true);
    std::cerr << smi << std::endl;
    TEST_ASSERT(smi == "CC[C@H](F)Cl");
    delete mol;
  }
  {
    std::string rdbase = getenv("RDBASE");
    std::string fName =
        rdbase + "/Code/GraphMol/FileParsers/test_data/basic_chiral1b.json";
    std::ifstream ifs(fName.c_str());

    RWMol *mol = JSONDataStreamToMol(ifs);
    TEST_ASSERT(mol);
    TEST_ASSERT(!mol->hasProp("_Name"));

    std::string smi = MolToSmiles(*mol, true);
    std::cerr << smi << std::endl;
    TEST_ASSERT(smi == "CC[C@H](F)Cl");
    delete mol;
  }
  BOOST_LOG(rdInfoLog) << "Finished" << std::endl;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  RDLog::InitLogs();

  std::string rdbase = getenv("RDBASE");

  testReadBasics();
  testWriteBasics();
  testRoundTripSmiles();
  testChiralBasics();

  return 0;
}
