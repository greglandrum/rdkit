//
//
//  Copyright (C) 2018-2019 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do
                           // this in one cpp file
#include "catch.hpp"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>

using namespace RDKit;

TEST_CASE("Github #1972", "[SMILES,bug]") {
  SECTION("basics") {
    std::vector<std::vector<std::string>> smiles = {
        {"[C@@]1(Cl)(F)(I).Br1", "[C@@](Br)(Cl)(F)(I)"},
        {"[C@@](Cl)(F)(I)1.Br1", "[C@@](Cl)(F)(I)Br"},
        {"[C@@](Cl)1(F)(I).Br1", "[C@@](Cl)(Br)(F)(I)"},
        {"[C@@](Cl)(F)1(I).Br1", "[C@@](Cl)(F)(Br)(I)"}};
    for (const auto &pr : smiles) {
      std::unique_ptr<ROMol> m1(SmilesToMol(pr[0]));
      std::unique_ptr<ROMol> m2(SmilesToMol(pr[1]));
      REQUIRE(m1);
      REQUIRE(m2);
      auto csmi1 = MolToSmiles(*m1);
      auto csmi2 = MolToSmiles(*m2);
      CHECK(csmi1 == csmi2);
    }
  }
  SECTION("further examples") {
    std::vector<std::vector<std::string>> smiles = {
        {"[C@@]1(Cl)2(I).Br1.F2", "[C@@](Br)(Cl)(F)(I)"},
        {"[C@@](Cl)2(I)1.Br1.F2", "[C@@](Cl)(F)(I)Br"},
        {"[C@@]12(Cl)(I).Br1.F2", "[C@@](Br)(F)(Cl)(I)"},
        {"[C@@]21(Cl)(I).Br1.F2", "[C@@](F)(Br)(Cl)(I)"},
        {"[C@@](Cl)12(I).Br1.F2", "[C@@](Cl)(Br)(F)(I)"},
        {"[C@@](Cl)21(I).Br1.F2", "[C@@](Cl)(F)(Br)(I)"},
        {"[C@@](Cl)(I)21.Br1.F2", "[C@@](Cl)(I)(F)(Br)"},
        {"[C@@](Cl)(I)12.Br1.F2", "[C@@](Cl)(I)(Br)(F)"}};
    for (const auto &pr : smiles) {
      std::unique_ptr<ROMol> m1(SmilesToMol(pr[0]));
      std::unique_ptr<ROMol> m2(SmilesToMol(pr[1]));
      REQUIRE(m1);
      REQUIRE(m2);
      auto csmi1 = MolToSmiles(*m1);
      auto csmi2 = MolToSmiles(*m2);
      CHECK(csmi1 == csmi2);
    }
  }
}

TEST_CASE("Github #2029", "[SMILES,bug]") {
  SECTION("wedging") {
    std::unique_ptr<ROMol> m1(SmilesToMol("CN[C@H](Cl)C(=O)O"));
    REQUIRE(m1);
    m1->getBondWithIdx(1)->setBondDir(Bond::BEGINWEDGE);
    bool doKekule = false, allBondsExplicit = false;
    CHECK("" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(1), -1, doKekule,
                                           allBondsExplicit));
    allBondsExplicit = true;
    CHECK("-" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(1), -1, doKekule,
                                            allBondsExplicit));
  }
  SECTION("direction") {
    std::unique_ptr<ROMol> m1(SmilesToMol("C/C=C/C"));
    REQUIRE(m1);
    bool doKekule = false, allBondsExplicit = false;
    CHECK("" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(0), -1, doKekule,
                                           allBondsExplicit));
    CHECK("" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(2), -1, doKekule,
                                           allBondsExplicit));
    allBondsExplicit = true;
    CHECK("/" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(0), -1, doKekule,
                                            allBondsExplicit));
    CHECK("/" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(2), -1, doKekule,
                                            allBondsExplicit));
  }
  SECTION("aromatic double bonds") {
    std::unique_ptr<RWMol> m1(SmilesToMol("c1ccccc1"));
    REQUIRE(m1);
    bool markAtomsBonds = false;
    MolOps::Kekulize(*m1, markAtomsBonds);
    bool doKekule = false, allBondsExplicit = false;
    CHECK("" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(0), -1, doKekule,
                                           allBondsExplicit));
    CHECK("" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(1), -1, doKekule,
                                           allBondsExplicit));
    allBondsExplicit = true;
    CHECK("=" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(0), -1, doKekule,
                                            allBondsExplicit));
    CHECK("-" == SmilesWrite::GetBondSmiles(m1->getBondWithIdx(1), -1, doKekule,
                                            allBondsExplicit));
  }
}

TEST_CASE("Smiles literals", "[SMILES]") {
  auto mol = "c1ccccc1"_smiles;
  REQUIRE(mol);
  CHECK(6 == mol->getNumAtoms());
  auto fail1 = "c1ccccc"_smiles;
  REQUIRE(!fail1);
  auto fail2 = "c1cccn1"_smiles;
  REQUIRE(!fail2);
}

TEST_CASE("Smarts literals", "[Smarts]") {
  auto mol = "c1ccc[c,n]c1"_smarts;
  REQUIRE(mol);
  CHECK(6 == mol->getNumAtoms());
  auto fail1 = "c1ccccc"_smarts;
  REQUIRE(!fail1);
  auto mol2 = "c1cccn1"_smarts;
  REQUIRE(mol2);
}

TEST_CASE(
    "github #2197 and #2237: handling of aromatic main group atoms in SMARTS",
    "[Smarts]") {
  std::vector<std::string> smarts = {
      "[si]1ccccc1",
      "[as]1ccccc1",
      "[se]1ccccc1",
      "[te]1ccccc1",

  };
  SECTION("#2197") {
    for (const auto sma : smarts) {
      std::unique_ptr<ROMol> mol(SmartsToMol(sma));
      REQUIRE(mol);
      CHECK(6 == mol->getNumAtoms());
      REQUIRE(mol->getAtomWithIdx(0)->hasQuery());
      REQUIRE(static_cast<QueryAtom *>(mol->getAtomWithIdx(0))
                  ->getQuery()
                  ->getDescription() == "AtomType");
    }
  }
  SECTION("#2237") {
    for (const auto sma : smarts) {
      std::unique_ptr<ROMol> mol(SmartsToMol(sma));
      REQUIRE(mol);
      REQUIRE(MolToSmarts(*mol) == sma);
    }
  }
}

TEST_CASE("Github #2148", "[bug, Smiles, Smarts]") {
  SECTION("SMILES") {
    auto mol = "C(=C\\F)\\4.O=C1C=4CCc2ccccc21"_smiles;
    REQUIRE(mol);
    REQUIRE(mol->getBondBetweenAtoms(0, 5));
    CHECK(mol->getBondBetweenAtoms(0, 5)->getBondType() == Bond::DOUBLE);
  }
  SECTION("SMILES edges") {
    auto m1 = "C/C=C/C"_smiles;
    REQUIRE(m1);
    CHECK(m1->getBondBetweenAtoms(2, 1)->getBondType() == Bond::DOUBLE);
    CHECK(m1->getBondBetweenAtoms(2, 1)->getStereo() != Bond::STEREONONE);

    {
      std::vector<std::string> smis = {"C1=C/C.C/1", "C/1=C/C.C1",
                                       "C-1=C/C.C/1", "C/1=C/C.C-1"};
      for (auto smi : smis) {
        std::unique_ptr<RWMol> mol(SmilesToMol(smi));
        REQUIRE(mol);
        CHECK(mol->getBondBetweenAtoms(0, 3)->getBondType() == Bond::SINGLE);
        CHECK(mol->getBondBetweenAtoms(0, 3)->getBondDir() != Bond::NONE);
        CHECK(mol->getBondBetweenAtoms(0, 1)->getBondType() == Bond::DOUBLE);
        CHECK(mol->getBondBetweenAtoms(0, 1)->getStereo() != Bond::STEREONONE);
      }
    }
  }

  SECTION("Writing SMILES") {
    auto mol = "C/C=c1/ncc(=C)cc1"_smiles;
    REQUIRE(mol);
    REQUIRE(mol->getBondBetweenAtoms(1, 2));
    CHECK(mol->getBondBetweenAtoms(1, 2)->getBondType() == Bond::DOUBLE);
    CHECK(mol->getBondBetweenAtoms(1, 2)->getStereo() == Bond::STEREOE);
    auto smi = MolToSmiles(*mol);
    CHECK(smi=="C=c1cc/c(=C\\C)nc1");    
  }


}
