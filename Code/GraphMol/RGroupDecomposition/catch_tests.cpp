//
//
//  Copyright (C) 2019 Greg Landrum and T5 Informatics GmbH
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
#include <GraphMol/RGroupDecomposition/RGroupDecomp.h>
#include <GraphMol/Depictor/RDDepictor.h>

using namespace RDKit;

TEST_CASE("A set of linked improvements", "[RGD]") {
  std::vector<std::string> smiles = {"c1c(F)cccn1", "c1c(Cl)c(C)ccn1",
                                     "c1c(F)c(C)ccn1"};
  std::vector<ROMOL_SPTR> mols;
  for (const auto &smi : smiles) {
    ROMOL_SPTR m(SmilesToMol(smi));
    REQUIRE(m);
    mols.push_back(m);
  }
  auto core = "c1c([*:2])c([*:1])ccn1"_smiles;
  REQUIRE(core);

  SECTION("github #2332: addHs() call should set coords") {
    // make sure we have coords on the molecules
    for (auto &mol : mols) {
      RDDepict::compute2DCoords(*mol);
    }
    RGroupDecompositionParameters params;
    params.removeHydrogensPostMatch = false;
    RGroupDecomposition rgd(*core, params);
    for (auto &mol : mols) {
      CHECK(rgd.add(*mol) >= 0);
    }
    CHECK(rgd.process());
    auto columns = rgd.getRGroupsAsColumns();
    for (auto col : columns) {
      if (col.first == "Core") continue;
      for (auto mol : col.second) {
        REQUIRE(mol);
        REQUIRE(mol->getNumConformers() == 1);
        const auto &conf = mol->getConformer();
        for (size_t i = 0; i < mol->getNumAtoms(); ++i) {
          // we don't really care what the coordinates are, just that something
          // has been set. :-)
          auto sum = abs(conf.getAtomPos(i).x) + abs(conf.getAtomPos(i).y);
          CHECK(sum != 0.0);
        }
      }
    }
  }
  SECTION("github #2331: default to removeHydrogensPostMatch=true") {
    RGroupDecomposition rgd(*core);
    for (auto &mol : mols) {
      CHECK(rgd.add(*mol) >= 0);
    }
    CHECK(rgd.process());
    auto columns = rgd.getRGroupsAsColumns();
    for (auto col : columns) {
      // if (col.first == "Core") continue;
      for (auto mol : col.second) {
        REQUIRE(mol);
        // skip r groups that are just Hs
        if (mol->getNumAtoms() == 2) continue;
        for (const auto atom : mol->atoms()) {
          CHECK(atom->getAtomicNum() != 1);
        }
      }
    }
  }
  SECTION(
      "github #2330: remove RGroups that are a mix of Hs and empty values") {
    auto core2 = "c1c([*:2])c([*:1])ncn1"_smiles;
    REQUIRE(core2);
    std::vector<ROMOL_SPTR> cores;
    cores.push_back(ROMOL_SPTR(new ROMol(*core)));
    cores.push_back(ROMOL_SPTR(new ROMol(*core2)));

    RGroupDecomposition rgd(cores);
    for (auto &mol : mols) {
      CHECK(rgd.add(*mol) >= 0);
    }
    // add an extra molecule
    auto m = "c1c(F)cncn1"_smiles;
    REQUIRE(m);
    std::cerr << "-----------------------------------" << std::endl;
    CHECK(rgd.add(*m) >= 0);
    CHECK(rgd.process());
    auto columns = rgd.getRGroupsAsColumns();
    for (const auto elem : columns) {
      std::cerr << elem.first << " " << elem.second.size() << std::endl;
    }
    CHECK(columns.size() == 3);
    CHECK(columns.find("Core") != columns.end());
    CHECK(columns.find("R1") != columns.end());
    CHECK(columns.find("R2") != columns.end());
    CHECK(columns.find("R3") == columns.end());
  }
}
