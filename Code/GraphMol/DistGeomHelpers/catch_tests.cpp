//
//  Copyright (C) 2021-2024 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/test.h>
#include <catch2/catch_all.hpp>

#include <RDGeneral/RDLog.h>
#include <GraphMol/test_fixtures.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Atropisomers.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ForceFieldHelpers/UFF/UFF.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/TorsionPreferences.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include "Embedder.h"
#include "BoundsMatrixBuilder.h"
#include <tuple>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>

using namespace RDKit;

TEST_CASE("Torsions not found in fused macrocycles", "[macrocycles]") {
  RDLog::InitLogs();
  SECTION("reported") {
    // this is 6VY8 from the PDB
    auto mol1 =
        "CC[C@H](C)[C@@H]1NC(=O)[C@@H]2CCCN2C(=O)[C@@H]2CCCN2C(=O)[C@H]([C@@H](C)CC)NC(=O)[C@H](CO)NC(=O)[C@H](CCCC[NH3+])NC(=O)[C@H]([C@@H](C)O)NC(O)[C@@H]2CN3NNC[C@H]3C[C@H](NC1=O)C(O)N[C@@H](Cc1ccccc1)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CC(=O)O)C(=O)NCC(=O)N[C@@H](CCCNC(N)=[NH2+])C(=O)N2"_smiles;
    REQUIRE(mol1);
    MolOps::addHs(*mol1);
    ForceFields::CrystalFF::CrystalFFDetails details;
    bool useExpTorsions = true;
    bool useSmallRingTorsions = false;
    bool useMacrocycleTorsions = true;
    bool useBasicKnowledge = true;
    unsigned int version = 2;
    bool verbose = true;
    std::stringstream sstrm;
    rdInfoLog->SetTee(sstrm);
    ForceFields::CrystalFF::getExperimentalTorsions(
        *mol1, details, useExpTorsions, useSmallRingTorsions,
        useMacrocycleTorsions, useBasicKnowledge, version, verbose);
    rdInfoLog->ClearTee();
    auto txt = sstrm.str();
    CHECK(txt.find("{9-}") != std::string::npos);
  }
  SECTION("edges") {
    std::vector<std::tuple<std::string, bool, unsigned int>> tests{
        {"O=C1CNC(=O)C2CCC(N1)NC(=O)CNC2=O", true, 15},  // 9-9
        {"O=C1NC2CCC(C(=O)N1)C(=O)NCC(=O)N2", true, 4},  // 9-8
        {"O=C1NC2CCC(C(=O)N1)C(=O)NC(=O)N2", false, 0},  // 8-8
        {"O=C1CC(=O)NC2NC(=O)CC(=O)NC(N1)NC(=O)CC(=O)N2", true,
         18}};  // 12-12-12
    for (const auto &tpl : tests) {
      std::unique_ptr<RWMol> m{SmilesToMol(std::get<0>(tpl))};
      REQUIRE(m);
      MolOps::addHs(*m);
      ForceFields::CrystalFF::CrystalFFDetails details;
      bool useExpTorsions = true;
      bool useSmallRingTorsions = false;
      bool useMacrocycleTorsions = true;
      bool useBasicKnowledge = true;
      unsigned int version = 2;
      bool verbose = true;
      std::stringstream sstrm;
      rdInfoLog->SetTee(sstrm);
      std::cerr << "-----------" << std::endl;
      ForceFields::CrystalFF::getExperimentalTorsions(
          *m, details, useExpTorsions, useSmallRingTorsions,
          useMacrocycleTorsions, useBasicKnowledge, version, verbose);
      rdInfoLog->ClearTee();
      auto txt = sstrm.str();
      if (std::get<1>(tpl)) {
        CHECK(txt.find("{9-}") != std::string::npos);
      } else {
        CHECK(txt.find("{9-}") == std::string::npos);
      }
      CHECK(details.expTorsionAngles.size() == std::get<2>(tpl));
    }
  }
}

namespace {
void compareConfs(const ROMol *m, const ROMol *expected, int molConfId = -1,
                  int expectedConfId = -1) {
  PRECONDITION(m, "bad pointer");
  PRECONDITION(expected, "bad pointer");
  TEST_ASSERT(m->getNumAtoms() == expected->getNumAtoms());
  const Conformer &conf1 = m->getConformer(molConfId);
  const Conformer &conf2 = expected->getConformer(expectedConfId);
  for (unsigned int i = 0; i < m->getNumAtoms(); i++) {
    TEST_ASSERT(m->getAtomWithIdx(i)->getAtomicNum() ==
                expected->getAtomWithIdx(i)->getAtomicNum());

    RDGeom::Point3D pt1i = conf1.getAtomPos(i);
    RDGeom::Point3D pt2i = conf2.getAtomPos(i);
    TEST_ASSERT((pt1i - pt2i).length() < 0.05);
  }
}
}  // namespace

TEST_CASE("update parameters from JSON") {
  std::string rdbase = getenv("RDBASE");
  SECTION("DG") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.dg.mol";
    std::unique_ptr<RWMol> ref{MolFileToMol(fname, true, false)};
    REQUIRE(ref);
    std::unique_ptr<RWMol> mol{SmilesToMol("OCCC")};
    REQUIRE(mol);
    MolOps::addHs(*mol);
    CHECK(ref->getNumAtoms() == mol->getNumAtoms());
    DGeomHelpers::EmbedParameters params;
    std::string json = R"JSON({"randomSeed":42})JSON";
    DGeomHelpers::updateEmbedParametersFromJSON(params, json);
    CHECK(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    compareConfs(ref.get(), mol.get());
  }
  SECTION("ETKDG") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.mol";
    std::unique_ptr<RWMol> ref{MolFileToMol(fname, true, false)};
    REQUIRE(ref);
    std::unique_ptr<RWMol> mol{SmilesToMol("OCCC")};
    REQUIRE(mol);
    MolOps::addHs(*mol);
    CHECK(ref->getNumAtoms() == mol->getNumAtoms());
    DGeomHelpers::EmbedParameters params;
    std::string json = R"JSON({"randomSeed":42,
    "useExpTorsionAnglePrefs":true,
    "useBasicKnowledge":true})JSON";
    DGeomHelpers::updateEmbedParametersFromJSON(params, json);
    CHECK(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    compareConfs(ref.get(), mol.get());
  }
  SECTION("ETKDGv2") {
    std::string fname =
        rdbase +
        "/Code/GraphMol/DistGeomHelpers/test_data/torsion.etkdg.v2.mol";
    std::unique_ptr<RWMol> ref{MolFileToMol(fname, true, false)};
    REQUIRE(ref);
    std::unique_ptr<RWMol> mol{SmilesToMol("n1cccc(C)c1ON")};
    REQUIRE(mol);
    MolOps::addHs(*mol);
    CHECK(ref->getNumAtoms() == mol->getNumAtoms());
    DGeomHelpers::EmbedParameters params;
    std::string json = R"JSON({"randomSeed":42,
    "useExpTorsionAnglePrefs":true,
    "useBasicKnowledge":true,
    "ETversion":2})JSON";
    DGeomHelpers::updateEmbedParametersFromJSON(params, json);
    CHECK(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
    compareConfs(ref.get(), mol.get());
  }

  SECTION("setting atommap") {
    std::unique_ptr<RWMol> mol{SmilesToMol("OCCC")};
    REQUIRE(mol);
    MolOps::addHs(*mol);
    {
      DGeomHelpers::EmbedParameters params;
      std::string json = R"JSON({"randomSeed":42,
    "coordMap":{"0":[0,0,0],"1":[0,0,1.5],"2":[0,1.5,1.5]}})JSON";
      DGeomHelpers::updateEmbedParametersFromJSON(params, json);
      CHECK(DGeomHelpers::EmbedMolecule(*mol, params) == 0);
      delete params.coordMap;
      auto conf = mol->getConformer();
      auto v1 = conf.getAtomPos(0) - conf.getAtomPos(1);
      auto v2 = conf.getAtomPos(2) - conf.getAtomPos(1);
      CHECK(v1.angleTo(v2) == Catch::Approx(M_PI / 2).margin(0.15));
    }
  }
}

TEST_CASE(
    "github #4346: Specified cis/trans stereo being ignored during "
    "conformation generation in macrocycles") {
  SECTION("basics 1") {
    auto m1 = "C1C/C=C/CCCCCCCC1"_smiles;
    REQUIRE(m1);
    CHECK(m1->getBondBetweenAtoms(2, 3)->getStereo() ==
          Bond::BondStereo::STEREOE);
    MolOps::addHs(*m1);
    DGeomHelpers::EmbedParameters params = DGeomHelpers::KDG;
    params.randomSeed = 0xf00d;
    CHECK(DGeomHelpers::EmbedMolecule(*m1, params) != -1);
    MolOps::assignStereochemistryFrom3D(*m1);
    CHECK(m1->getBondBetweenAtoms(2, 3)->getStereo() ==
          Bond::BondStereo::STEREOE);
  }
  SECTION("basics 2") {
    auto m1 = "C1C/C=C\\CCCCCCCC1"_smiles;
    REQUIRE(m1);
    CHECK(m1->getBondBetweenAtoms(2, 3)->getStereo() ==
          Bond::BondStereo::STEREOZ);
    MolOps::addHs(*m1);
    DGeomHelpers::EmbedParameters params = DGeomHelpers::KDG;
    params.randomSeed = 0xf00d;
    CHECK(DGeomHelpers::EmbedMolecule(*m1, params) != -1);
    MolOps::assignStereochemistryFrom3D(*m1);
    CHECK(m1->getBondBetweenAtoms(2, 3)->getStereo() ==
          Bond::BondStereo::STEREOZ);
  }
}
TEST_CASE("nontetrahedral stereo", "[nontetrahedral]") {
  SECTION("bounds matrix basics") {
    {
      auto m = "Cl[Pt@SP1]([35Cl])([36Cl])[37Cl]"_smiles;
      REQUIRE(m);
      CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                           m->getAtomWithIdx(0))
                ->getIdx() == 3);
      CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                           m->getAtomWithIdx(2))
                ->getIdx() == 4);
      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(3)),
          Catch::Matchers::WithinAbs(180, 0.001));

      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(2)),
          Catch::Matchers::WithinAbs(90, 0.001));

      DistGeom::BoundsMatPtr bm{new DistGeom::BoundsMatrix(m->getNumAtoms())};
      DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
      DGeomHelpers::setTopolBounds(*m, bm);
      // std::cerr << *bm << std::endl;
      CHECK(bm->getLowerBound(0, 3) - bm->getLowerBound(0, 2) > 1.0);
      CHECK(bm->getUpperBound(0, 3) - bm->getUpperBound(0, 2) > 1.0);
    }

    {
      auto m = "Cl[Pt@SP1]([35Cl])[36Cl]"_smiles;
      REQUIRE(m);
      CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                           m->getAtomWithIdx(0))
                ->getIdx() == 3);
      CHECK(!Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                            m->getAtomWithIdx(2)));
      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(3)),
          Catch::Matchers::WithinAbs(180, 0.001));

      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(2)),
          Catch::Matchers::WithinAbs(90, 0.001));

      DistGeom::BoundsMatPtr bm{new DistGeom::BoundsMatrix(m->getNumAtoms())};
      DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
      DGeomHelpers::setTopolBounds(*m, bm);
      // std::cerr << *bm << std::endl;
      CHECK(bm->getLowerBound(0, 3) - bm->getLowerBound(0, 2) > 1.0);
      CHECK(bm->getUpperBound(0, 3) - bm->getUpperBound(0, 2) > 1.0);
    }

    {
      // note that things aren't quite as nice here since we don't actually have
      // TBP UFF parameters
      auto m = "Cl[Pt@TB1]([35Cl])([36Cl])([37Cl])[38Cl]"_smiles;
      REQUIRE(m);
      CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                           m->getAtomWithIdx(0))
                ->getIdx() == 5);
      CHECK(!Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                            m->getAtomWithIdx(2)));
      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(5)),
          Catch::Matchers::WithinAbs(180, 0.001));

      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(2)),
          Catch::Matchers::WithinAbs(90, 0.001));
      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(3), m->getAtomWithIdx(2)),
          Catch::Matchers::WithinAbs(120, 0.001));

      DistGeom::BoundsMatPtr bm{new DistGeom::BoundsMatrix(m->getNumAtoms())};
      DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
      DGeomHelpers::setTopolBounds(*m, bm);
      CHECK(bm->getLowerBound(0, 5) - bm->getLowerBound(0, 2) > 0.5);
      CHECK(bm->getUpperBound(0, 5) - bm->getUpperBound(0, 2) > 0.5);
      CHECK(bm->getLowerBound(0, 5) - bm->getLowerBound(2, 3) > 0.5);
      CHECK(bm->getUpperBound(0, 5) - bm->getUpperBound(2, 3) > 0.5);
      CHECK(bm->getLowerBound(2, 3) - bm->getLowerBound(0, 2) > 0.5);
      CHECK(bm->getUpperBound(2, 3) - bm->getUpperBound(0, 2) > 0.5);
    }
    {
      auto m = "Cl[Th@OH1]([35Cl])([36Cl])([37Cl])([38Cl])[39Cl]"_smiles;
      REQUIRE(m);
      CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                           m->getAtomWithIdx(0))
                ->getIdx() == 6);
      CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                           m->getAtomWithIdx(2))
                ->getIdx() == 4);
      CHECK(Chirality::getChiralAcrossAtom(m->getAtomWithIdx(1),
                                           m->getAtomWithIdx(3))
                ->getIdx() == 5);

      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(6)),
          Catch::Matchers::WithinAbs(180, 0.001));

      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(0), m->getAtomWithIdx(2)),
          Catch::Matchers::WithinAbs(90, 0.001));
      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(4), m->getAtomWithIdx(2)),
          Catch::Matchers::WithinAbs(180, 0.001));
      CHECK_THAT(
          Chirality::getIdealAngleBetweenLigands(
              m->getAtomWithIdx(1), m->getAtomWithIdx(3), m->getAtomWithIdx(2)),
          Catch::Matchers::WithinAbs(90, 0.001));

      DistGeom::BoundsMatPtr bm{new DistGeom::BoundsMatrix(m->getNumAtoms())};
      DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
      DGeomHelpers::setTopolBounds(*m, bm);
      CHECK(bm->getLowerBound(0, 6) - bm->getLowerBound(0, 2) > 0.5);
      CHECK(bm->getUpperBound(0, 6) - bm->getUpperBound(0, 3) > 0.5);
      CHECK(bm->getLowerBound(0, 6) - bm->getLowerBound(2, 3) > 0.5);
      CHECK(bm->getUpperBound(0, 6) - bm->getUpperBound(2, 4) < 0.01);
      CHECK(bm->getLowerBound(2, 4) - bm->getLowerBound(2, 3) > 0.5);
    }
  }
#if 1
  SECTION("Embedding") {
    {
      auto m = "Cl[Pt@SP1](<-N)(<-N)[Cl]"_smiles;
      REQUIRE(m);
      m->setProp("_Name", "cis platin");
      MolOps::addHs(*m);
      CHECK(DGeomHelpers::EmbedMolecule(*m) == 0);
      auto mb = MolToV3KMolBlock(*m);
      // std::cerr << mb << std::endl;
      std::unique_ptr<RWMol> m2(MolBlockToMol(mb));
      MolOps::assignStereochemistryFrom3D(*m2);
      CHECK(m2->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_SQUAREPLANAR);
      unsigned int perm = 100;
      CHECK(m2->getAtomWithIdx(1)->getPropIfPresent(
          common_properties::_chiralPermutation, perm));
      CHECK(perm == 1);
    }
    {
      auto m = "Cl[Pt@SP3](<-N)(<-N)[Cl]"_smiles;
      REQUIRE(m);
      m->setProp("_Name", "trans platin");
      MolOps::addHs(*m);
      CHECK(DGeomHelpers::EmbedMolecule(*m) == 0);
      auto mb = MolToV3KMolBlock(*m);
      // std::cerr << mb << std::endl;
      std::unique_ptr<RWMol> m2(MolBlockToMol(mb));
      MolOps::assignStereochemistryFrom3D(*m2);
      CHECK(m2->getAtomWithIdx(1)->getChiralTag() ==
            Atom::ChiralType::CHI_SQUAREPLANAR);
      unsigned int perm = 100;
      CHECK(m2->getAtomWithIdx(1)->getPropIfPresent(
          common_properties::_chiralPermutation, perm));
      CHECK(perm == 3);
    }
  }
#endif
}

TEST_CASE("problems with bounds matrix smoothing and aromatic sulfur") {
  SECTION("basics") {
    auto core = R"CTAB(test structure - renumbered
     RDKit          3D

  7  7  0  0  0  0  0  0  0  0999 V2000
   48.6842  -14.8137    0.1450 C   0  0  0  0  0  0  0  0  0  0  0  0
   48.0829  -13.5569    0.6868 C   0  0  0  0  0  0  0  0  0  0  0  0
   48.0162  -12.0909   -0.1327 S   0  0  0  0  0  0  0  0  0  0  0  0
   47.1565  -11.3203    1.0899 C   0  0  0  0  0  0  0  0  0  0  0  0
   46.9350  -12.2470    2.1088 C   0  0  0  0  0  0  0  0  0  0  0  0
   46.1942  -11.9293    3.3432 C   0  0  0  0  0  0  0  0  0  0  0  0
   47.4440  -13.4879    1.8745 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  7  2  0
  2  3  1  0
  7  5  1  0
  5  4  2  0
  5  6  1  0
  4  3  1  0
M  END)CTAB"_ctab;
    REQUIRE(core);
    auto thiaz = "Cc1scc(C)n1"_smiles;
    REQUIRE(thiaz);
    MolOps::addHs(*thiaz);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    const auto conf = core->getConformer();
    std::map<int, RDGeom::Point3D> cmap;
    for (unsigned i = 0; i < core->getNumAtoms(); ++i) {
      cmap[i] = conf.getAtomPos(i);
    }
    ps.coordMap = &cmap;
    ps.randomSeed = 0xf00d;
    auto cid = DGeomHelpers::EmbedMolecule(*thiaz, ps);
    CHECK(cid >= 0);
  }
  SECTION("bulk") {
    // run a bunch of molecules with S-containing aromatic heterocycles
    std::vector<std::string> smileses = {
        "[O-][S+](c1ccccn1)c1cncs1",
        "Cn1cccc1C(=O)Nc1nccs1",
        "Cc1csc(=N)n1-c1ccc(Cl)cc1",
        "Nc1ncc([S+]([O-])c2ncccn2)s1",
        "CCCN1CCC=C(c2csc(N)n2)C1",
        "CNc1ncc([S+]([O-])c2ccccn2)s1",
        "Cn1nnnc1SCc1nc2ccccc2s1",
        "CCCC(C(=O)Nc1nccs1)c1ccccc1",
        "Cc1ccc(NC(=O)c2sc(Cl)nc2C)c(C)c1",
        "CCc1nc(-c2ccc(Cl)cc2)sc1C(=O)OC",
        "Cc1nc(CNS(=O)(=O)c2ccc(Cl)cc2)cs1",
        "Cc1ccc2sc(C)[n+](CCC(C)S(=O)(=O)[O-])c2c1",
        "Nc1nc2c(s1)-c1ccccc1Sc1ccccc1-2",
        "COc1ccccc1OCC(=O)Nc1nc(C)c(C)s1",
        "COc1ccc(NC(=O)Nc2sc(=S)n(C)c2C)cc1",
        "C=CCNc1nc(-c2c[nH]c3c(CC)cccc23)cs1",
    };
    auto patt = "[s]1*c[!#6]c1"_smarts;
    REQUIRE(patt);
    for (const auto &smi : smileses) {
      INFO(smi);
      std::unique_ptr<RWMol> mol{SmilesToMol(smi)};
      REQUIRE(mol);
      MolOps::addHs(*mol);
      DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
      ps.randomSeed = 0xf00d;
      auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
      REQUIRE(cid >= 0);
      UFF::UFFOptimizeMolecule(*mol);

      auto match = SubstructMatch(*mol, *patt);
      REQUIRE(match.size() >= 1);

      const auto conf = mol->getConformer();
      std::map<int, RDGeom::Point3D> cmap;
      for (auto &mi : match[0]) {
        cmap[mi.second] = conf.getAtomPos(mi.second);
      }
      ps.coordMap = &cmap;
      auto cid2 = DGeomHelpers::EmbedMolecule(*mol, ps);
      CHECK(cid2 >= 0);
    }
  }
  SECTION("phosphorous") {
    std::vector<std::string> smileses = {
        "CCOC(=O)c1pc(P(Cl)Cl)c2n1[C@@H](C)C(=O)Nc1ccc(C)cc1-2",
        "N(c1c(O)ccc2c(P(Cl)Cl)pc(C(=O)O)n12)[N+](=O)[O-]",
    };
    auto patt = "[p]1*c[!#6]c1"_smarts;
    REQUIRE(patt);
    for (const auto &smi : smileses) {
      INFO(smi);
      std::unique_ptr<RWMol> mol{SmilesToMol(smi)};
      REQUIRE(mol);
      MolOps::addHs(*mol);
      DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
      ps.randomSeed = 0xf00d;
      auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
      REQUIRE(cid >= 0);
      UFF::UFFOptimizeMolecule(*mol);

      auto match = SubstructMatch(*mol, *patt);
      REQUIRE(match.size() >= 1);

      const auto conf = mol->getConformer();
      std::map<int, RDGeom::Point3D> cmap;
      for (auto &mi : match[0]) {
        cmap[mi.second] = conf.getAtomPos(mi.second);
      }
      ps.coordMap = &cmap;
      auto cid2 = DGeomHelpers::EmbedMolecule(*mol, ps);
      CHECK(cid2 >= 0);
    }
  }
}

TEST_CASE("double bond stereo not honored in conformer generator") {
  SECTION("mol 1 basics") {
    // this test used to fail
    // from the platinum set
    auto m = "O=C1OCC/C=C/CC/C=C/C(=N/OCC(=O)N2CCCCC2)Cc2cc(O)cc(O)c21"_smiles;
    REQUIRE(m);
    RWMol cp(*m);
    MolOps::addHs(cp);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.randomSeed = 0xf00d + 81;
    auto cid = DGeomHelpers::EmbedMolecule(cp, ps);
    REQUIRE(cid >= 0);
    MolOps::assignStereochemistryFrom3D(cp);
    // std::cerr << MolToMolBlock(cp) << std::endl;
    for (const auto bnd : cp.bonds()) {
      if (bnd->getBondType() == Bond::BondType::DOUBLE) {
        INFO(bnd->getIdx());
        CHECK(bnd->getStereo() ==
              m->getBondWithIdx(bnd->getIdx())->getStereo());
      }
    }
  }
  SECTION("mol 1 multiple loops") {
    // from the platinum set
    auto m = "O=C1OCC/C=C/CC/C=C/C(=N/OCC(=O)N2CCCCC2)Cc2cc(O)cc(O)c21"_smiles;
    REQUIRE(m);
    RWMol cp(*m);
    MolOps::addHs(cp);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    for (unsigned int iter = 0; iter < 10; ++iter) {
      RWMol lcp(cp);
      ps.randomSeed = 0xf00d + iter;
      auto cid = DGeomHelpers::EmbedMolecule(lcp, ps);
      REQUIRE(cid >= 0);
      MolOps::assignStereochemistryFrom3D(lcp);
      // std::cerr << MolToMolBlock(cp) << std::endl;
      for (const auto bnd : lcp.bonds()) {
        if (bnd->getBondType() == Bond::BondType::DOUBLE) {
          INFO(iter);
          CHECK(bnd->getStereo() ==
                m->getBondWithIdx(bnd->getIdx())->getStereo());
        }
      }
    }
  }

  SECTION("github #5913") {
    auto m = "[H]/C(F)=C([H])\\C([H])=C(/[H])Br"_smiles;
    REQUIRE(m);

    RWMol cp(*m);
    MolOps::addHs(cp);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    for (unsigned int iter = 0; iter < 50; ++iter) {
      RWMol lcp(cp);
      ps.randomSeed = 0 + iter;
      auto cid = DGeomHelpers::EmbedMolecule(lcp, ps);
      REQUIRE(cid >= 0);
      MolOps::assignStereochemistryFrom3D(lcp);
      // std::cerr << MolToMolBlock(cp) << std::endl;
      for (const auto bnd : lcp.bonds()) {
        if (bnd->getBondType() == Bond::BondType::DOUBLE) {
          INFO(iter);
          CHECK(bnd->getStereo() ==
                m->getBondWithIdx(bnd->getIdx())->getStereo());
        }
      }
    }
  }

  SECTION("github #5283") {
    UseLegacyStereoPerceptionFixture useLegacy(false);
    auto m =
        "Cc3nn(CC(=O)N2CCN(c1ccccc1)CC2)c(C)c3/N=N\\c6ccc(CNC(=O)CCC(=O)Nc4cccc5C(=O)NCc45)cc6"_smiles;
    REQUIRE(m);
    RWMol cp(*m);
    MolOps::addHs(cp);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.enforceChirality = false;
    for (unsigned int iter = 0; iter < 10; ++iter) {
      INFO(iter);
      RWMol lcp(cp);
      ps.randomSeed = 140 + iter;
      auto cid = DGeomHelpers::EmbedMolecule(lcp, ps);
      REQUIRE(cid >= 0);
      MolOps::assignStereochemistryFrom3D(lcp, cid, true);
      auto bnd = lcp.getBondBetweenAtoms(22, 23);
      REQUIRE(bnd);
      REQUIRE(bnd->getBondType() == Bond::BondType::DOUBLE);
      CHECK(bnd->getStereo() == m->getBondWithIdx(bnd->getIdx())->getStereo());
    }
  }
}

TEST_CASE("tracking failure causes"){SECTION("basics"){
    auto mol =
        "C=CC1=C(N)Oc2cc1c(-c1cc(C(C)O)cc(=O)cc1C1NCC(=O)N1)c(OC)c2OC"_smiles;
REQUIRE(mol);
MolOps::addHs(*mol);
DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
ps.randomSeed = 0xf00d;
ps.trackFailures = true;
ps.maxIterations = 50;
ps.randomSeed = 42;
auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
CHECK(cid < 0);

CHECK(ps.failures[DGeomHelpers::EmbedFailureCauses::INITIAL_COORDS] > 5);
CHECK(ps.failures[DGeomHelpers::EmbedFailureCauses::ETK_MINIMIZATION] > 10);

auto fail_cp = ps.failures;
// make sure we reset the counts each time
cid = DGeomHelpers::EmbedMolecule(*mol, ps);
CHECK(ps.failures == fail_cp);
}
SECTION("chirality") {
  auto mol = R"CTAB(
  Ketcher  1102315302D 1   1.00000     0.00000     0

 10 11  0  0  1  0  0  0  0  0999 V2000
   10.1340  -11.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.1340  -12.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.0000  -12.5250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.8660  -12.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.8660  -11.0250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.0000  -10.5250    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   11.0000  -11.5250    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   11.2588  -12.4909    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.2680  -10.5250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.7629  -12.4673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  6  1  0     0  0
  1  2  1  0     0  0
  2  3  1  0     0  0
  3  4  1  0     0  0
  4  5  1  0     0  0
  5  6  1  0     0  0
  1  7  1  0     0  0
  7  8  1  0     0  0
  8  4  1  0     0  0
  1  9  1  1     0  0
  4 10  1  1     0  0
M  END
)CTAB"_ctab;
  REQUIRE(mol);
  MolOps::addHs(*mol);
  DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
  ps.randomSeed = 0xf00d;
  ps.trackFailures = true;
  ps.maxIterations = 50;
  auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
  CHECK(cid < 0);
  CHECK(ps.failures[DGeomHelpers::EmbedFailureCauses::INITIAL_COORDS] > 5);
  CHECK(ps.failures[DGeomHelpers::EmbedFailureCauses::FINAL_CHIRAL_BOUNDS] > 5);
}

#ifdef RDK_TEST_MULTITHREADED
SECTION("multithreaded") {
  auto mol =
      "C=CC1=C(N)Oc2cc1c(-c1cc(C(C)O)cc(=O)cc1C1NCC(=O)N1)c(OC)c2OC"_smiles;
  REQUIRE(mol);
  MolOps::addHs(*mol);
  DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
  ps.randomSeed = 0xf00d;
  ps.trackFailures = true;
  ps.maxIterations = 10;
  ps.randomSeed = 42;
  auto cids = DGeomHelpers::EmbedMultipleConfs(*mol, 20, ps);

  DGeomHelpers::EmbedParameters ps2 = ps;
  ps2.numThreads = 4;

  auto cids2 = DGeomHelpers::EmbedMultipleConfs(*mol, 20, ps2);
  CHECK(cids2 == cids);

  CHECK(ps.failures == ps2.failures);
}
#endif
}

TEST_CASE("Github #5883: confgen failing for chiral N in a three ring") {
  SECTION("basics1") {
    auto mol = "N1[C@H-]C1"_smiles;
    REQUIRE(mol);
    MolOps::addHs(*mol);
    mol->getAtomWithIdx(1)->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.randomSeed = 42;
    ps.maxIterations = 1;
    auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
    CHECK(cid >= 0);
  }
  SECTION("basics2") {
    auto mol = "N1[N@H]C1"_smiles;
    REQUIRE(mol);
    MolOps::addHs(*mol);
    mol->getAtomWithIdx(1)->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.randomSeed = 42;
    ps.maxIterations = 1;
    auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
    CHECK(cid >= 0);
  }
  SECTION("no ring") {
    auto mol = "N[C@H-]C"_smiles;
    REQUIRE(mol);
    MolOps::addHs(*mol);
    mol->getAtomWithIdx(1)->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.randomSeed = 42;
    ps.maxIterations = 1;
    auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
    CHECK(cid >= 0);
  }
}

TEST_CASE("Github #6365: cannot generate conformers for PF6- or SF6") {
  SECTION("basics") {
    std::vector<std::string> smileses = {"S(F)(F)(F)(F)(F)F",
                                         "[P-](F)(F)(F)(F)(F)F"};
    for (const auto &smi : smileses) {
      std::unique_ptr<RWMol> mol{SmilesToMol(smi)};
      REQUIRE(mol);
      DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
      ps.randomSeed = 42;
      ps.useRandomCoords = true;
      auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
      CHECK(cid >= 0);
    }
  }
}

TEST_CASE("Sequential random seeds") {
  SECTION("basics") {
    auto mol = "CCCCCCCCCCCC"_smiles;
    REQUIRE(mol);
    MolOps::addHs(*mol);

    RWMol mol2(*mol);

    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.enableSequentialRandomSeeds = true;
    ps.useRandomCoords = true;
    ps.randomSeed = 0xf00d;
    auto cids = DGeomHelpers::EmbedMultipleConfs(*mol, 10, ps);
    CHECK(cids.size() == 10);
    ps.randomSeed = 0xf00d + 5;
    auto cids2 = DGeomHelpers::EmbedMultipleConfs(mol2, 5, ps);
    CHECK(cids2.size() == 5);

    compareConfs(mol.get(), &mol2, 5, 0);
  }
}

TEST_CASE("Macrocycle bounds matrix") {
  SECTION("basics") {
    auto mol = "C1/C=C/C=C/CCCCCCCCC1"_smiles;
    REQUIRE(mol);
    MolOps::addHs(*mol);

    DistGeom::BoundsMatPtr bm{new DistGeom::BoundsMatrix(mol->getNumAtoms())};
    DGeomHelpers::initBoundsMat(bm, 0.0, 1000.0);
    DGeomHelpers::setTopolBounds(*mol, bm, true, false, true);
    CHECK(bm->getLowerBound(1, 18) > 2.6);
    CHECK(bm->getLowerBound(1, 18) < 2.7);
    CHECK(bm->getLowerBound(4, 17) > 2.6);
    CHECK(bm->getLowerBound(4, 17) < 2.7);

    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.randomSeed = 0;

    auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
    CHECK(cid >= 0);
    const auto conf = mol->getConformer(cid);
    RDGeom::Point3D pos_1 = conf.getAtomPos(1);
    RDGeom::Point3D pos_4 = conf.getAtomPos(4);
    CHECK((pos_1 - pos_4).length() < 3.6);
    CHECK((pos_1 - pos_4).length() > 3.5);
  }
}

TEST_CASE("atropisomers and embedding") {
  SECTION("basics") {
    auto mol =
        "Cc1cccc(O)c1-c1c(N)cccc1Cl |(-8.88571,2.09707,;-8.17143,3.33425,;-6.74286,3.33425,;-6.02857,4.57143,;-6.74286,5.80861,;-8.17143,5.80861,;-8.88571,7.04579,;-8.88571,4.57143,;-10.3143,4.57143,;-11.0286,5.80861,;-10.3143,7.04579,;-12.4571,5.80861,;-13.1714,4.57143,;-12.4571,3.33425,;-11.0286,3.33425,;-10.3143,2.09707,),wU:8.15|"_smiles;
    REQUIRE(mol);
    REQUIRE(mol->getBondWithIdx(7)->getBondType() == Bond::BondType::SINGLE);
    REQUIRE(mol->getBondWithIdx(7)->getStereo() ==
            Bond::BondStereo::STEREOATROPCCW);
    MolOps::addHs(*mol);
    // mol->debugMol(std::cerr);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.randomSeed = 0xf00d;
    {
      auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
      REQUIRE(cid >= 0);
      const auto conf = mol->getConformer(cid);

      Atropisomers::AtropAtomAndBondVec abvs[2];
      REQUIRE(Atropisomers::getAtropisomerAtomsAndBonds(mol->getBondWithIdx(7),
                                                        abvs, *mol));
      auto pos_1 = conf.getAtomPos(7);
      auto pos_2 = conf.getAtomPos(8);
      auto pos_3 = conf.getAtomPos(1);
      auto pos_4 = conf.getAtomPos(9);
      auto v2 = pos_2 - pos_1;
      auto v3 = pos_3 - pos_1;
      auto v4 = pos_4 - pos_1;
      auto chiralVol = v3.crossProduct(v4).dotProduct(v2);
      CHECK(chiralVol < 0);
    }
    {
      RWMol mol2(*mol);
      mol2.getBondWithIdx(7)->setStereo(Bond::BondStereo::STEREOATROPCW);

      auto cid = DGeomHelpers::EmbedMolecule(mol2, ps);
      REQUIRE(cid >= 0);
      const auto conf = mol2.getConformer(cid);

      Atropisomers::AtropAtomAndBondVec abvs[2];
      REQUIRE(Atropisomers::getAtropisomerAtomsAndBonds(mol2.getBondWithIdx(7),
                                                        abvs, mol2));
      auto pos_1 = conf.getAtomPos(7);
      auto pos_2 = conf.getAtomPos(8);
      auto pos_3 = conf.getAtomPos(1);
      auto pos_4 = conf.getAtomPos(9);
      auto v2 = pos_2 - pos_1;
      auto v3 = pos_3 - pos_1;
      auto v4 = pos_4 - pos_1;
      auto chiralVol = v3.crossProduct(v4).dotProduct(v2);
      CHECK(chiralVol > 0);
    }
  }
}

TEST_CASE("atropisomers bulk") {
  std::string rdbase = getenv("RDBASE");
  std::string fname =
      rdbase + "/Code/GraphMol/DistGeomHelpers/test_data/atropisomers.sdf";
  SDMolSupplier sdsup(fname);

  auto params = DGeomHelpers::ETKDGv3;
  params.randomSeed = 0xf00d + 1;

  for (auto i = 0u; i < sdsup.length(); ++i) {
    std::unique_ptr<RWMol> mol(static_cast<RWMol *>(sdsup[i]));
    REQUIRE(mol);
    auto bondIdx = mol->getProp<unsigned int>("atrop bond");
    REQUIRE((mol->getBondWithIdx(bondIdx)->getStereo() ==
                 Bond::BondStereo::STEREOATROPCCW ||
             mol->getBondWithIdx(bondIdx)->getStereo() ==
                 Bond::BondStereo::STEREOATROPCW));
    auto atropInfo = mol->getProp<std::string>("atrop volume");
    std::vector<std::string> tokens;
    boost::split(tokens, atropInfo, boost::is_any_of(" \t"));
    REQUIRE(tokens.size() == 5);
    std::vector<unsigned int> atropAtoms(4);
    for (auto j = 0u; j < 4u; ++j) {
      atropAtoms[j] = std::stol(tokens[j]);
    }
    int vol = std::stol(tokens[4]);

    MolOps::addHs(*mol);
    unsigned int nconfs = 20;
    {
      auto cids = DGeomHelpers::EmbedMultipleConfs(*mol, nconfs, params);
      CHECK(cids.size() == nconfs);
      for (auto cid : cids) {
        const auto conf = mol->getConformer(cid);
        std::vector<RDGeom::Point3D> pts;
        for (auto idx : atropAtoms) {
          pts.push_back(conf.getAtomPos(idx));
        }
        auto v2 = pts[1] - pts[0];
        auto v3 = pts[2] - pts[0];
        auto v4 = pts[3] - pts[0];
        auto chiralVol = v3.crossProduct(v4).dotProduct(v2);
        INFO(cid << MolToV3KMolBlock(*mol, true, cid));
        CHECK(chiralVol * vol > 0);
      }
    }  // now swap the stereo and see if it still works
    mol->getBondWithIdx(bondIdx)->setStereo(
        mol->getBondWithIdx(bondIdx)->getStereo() ==
                Bond::BondStereo::STEREOATROPCCW
            ? Bond::BondStereo::STEREOATROPCW
            : Bond::BondStereo::STEREOATROPCCW);
    {
      auto cids = DGeomHelpers::EmbedMultipleConfs(*mol, nconfs, params);
      CHECK(cids.size() == nconfs);
      for (auto cid : cids) {
        const auto conf = mol->getConformer(cid);
        std::vector<RDGeom::Point3D> pts;
        for (auto idx : atropAtoms) {
          pts.push_back(conf.getAtomPos(idx));
        }
        auto v2 = pts[1] - pts[0];
        auto v3 = pts[2] - pts[0];
        auto v4 = pts[3] - pts[0];
        auto chiralVol = v3.crossProduct(v4).dotProduct(v2);
        INFO(cid << MolToV3KMolBlock(*mol, true, cid));
        CHECK(chiralVol * vol < 0);
      }
    }
  }
}

TEST_CASE(
    "Github #7109: wrong stereochemistry in ring from stereospecific SMILES") {
  SECTION("basics") {
    auto m = "C1[C@H](C#CC#C)CC[C@H](C#CC#C)C1"_smiles;
    REQUIRE(m);
    MolOps::addHs(*m);
    REQUIRE(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
    REQUIRE(m->getAtomWithIdx(8)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);

    DGeomHelpers::EmbedParameters ps = DGeomHelpers::KDG;
    {  // this always worked
      ps.randomSeed = 0xC0FFEE;
      auto cid = DGeomHelpers::EmbedMolecule(*m, ps);
      CHECK(cid >= 0);
      MolOps::assignStereochemistryFrom3D(*m, cid);
      CHECK(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
      CHECK(m->getAtomWithIdx(8)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
    }
    {  // this failed
      ps.randomSeed = 0xC0FFEE + 123;
      auto cid = DGeomHelpers::EmbedMolecule(*m, ps);
      CHECK(cid >= 0);
      MolOps::assignStereochemistryFrom3D(*m, cid);
      CHECK(m->getAtomWithIdx(1)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
      CHECK(m->getAtomWithIdx(8)->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW);
    }
  }
}

TEST_CASE("user-provided torsions") {
  SECTION("basics") {
    auto mol = "O=CC=N"_smiles;
    REQUIRE(mol);
    // MolOps::addHs(*mol);
    DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
    ps.randomSeed = 0xf00d;
    ps.useRandomCoords = true;
    ps.useExpTorsionAnglePrefs = true;
    auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
    CHECK(cid >= 0);
    CHECK_THAT(MolTransforms::getDihedralDeg(mol->getConformer(), 0, 1, 2, 3),
               Catch::Matchers::WithinAbs(0.0, 0.1));
    // std::cerr << MolToV3KMolBlock(*mol) << std::endl;

    {  // silly test to make sure we can set the torsion
      DGeomHelpers::UserProvidedTorsion ut = {
          {0, 1, 2, 3}, {0.0, 15.0, 15.0, 0.0, 0.0, 0.0}, {1, -1, 1, 1, 1, 1}};
      ps.explicitTorsions.push_back(ut);
      cid = DGeomHelpers::EmbedMolecule(*mol, ps);
      CHECK(cid >= 0);
      CHECK_THAT(MolTransforms::getDihedralDeg(mol->getConformer(), 0, 1, 2, 3),
                 Catch::Matchers::WithinAbs(-46, 0.1));
    }

    {  // set it back to what ETKDG assigns
      DGeomHelpers::UserProvidedTorsion ut = {
          {0, 1, 2, 3}, {0.0, 15.0, 0.0, 0.0, 0.0, 0.0}, {1, -1, 1, 1, 1, 1}};
      ps.explicitTorsions.clear();
      ps.explicitTorsions.push_back(ut);
      cid = DGeomHelpers::EmbedMolecule(*mol, ps);
      CHECK(cid >= 0);
      CHECK_THAT(MolTransforms::getDihedralDeg(mol->getConformer(), 0, 1, 2, 3),
                 Catch::Matchers::WithinAbs(0, 0.1));
    }
    // std::cerr << MolToV3KMolBlock(*mol) << std::endl;
  }
}

TEST_CASE("WOT") {
  auto ctab = R"CTAB(8101969
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 55 58 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O 6.112500 -8.335700 -9.256800 0
M  V30 2 N 4.815500 -8.444700 -8.590700 0
M  V30 3 C 6.152100 -9.378400 -10.249000 0 CFG=2
M  V30 4 C 3.773400 -7.905600 -9.486900 0 CFG=1
M  V30 5 C 4.993700 -7.617200 -7.395700 0
M  V30 6 C 7.523000 -9.360100 -10.877200 0
M  V30 7 C 5.076100 -9.090000 -11.264900 0
M  V30 8 C 4.044400 -8.333500 -10.908800 0
M  V30 9 C 2.398300 -8.489100 -9.097700 0 CFG=2
M  V30 10 C 5.902400 -8.278800 -6.400000 0
M  V30 11 C 8.281000 -10.526100 -10.901200 0
M  V30 12 C 8.026700 -8.208700 -11.480300 0
M  V30 13 O 3.067800 -7.844500 -11.739700 0
M  V30 14 O 1.984000 -8.023900 -7.812300 0
M  V30 15 C 2.377000 -10.034600 -8.963800 0
M  V30 16 C 6.892000 -7.555500 -5.748200 0
M  V30 17 C 5.735400 -9.620900 -6.073300 0
M  V30 18 C 9.528400 -10.549400 -11.524700 0
M  V30 19 C 9.279300 -8.227200 -12.080500 0
M  V30 20 C 3.272300 -8.077000 -13.136600 0
M  V30 21 C 1.198900 -9.079200 -7.243200 0
M  V30 22 O 1.869100 -10.273000 -7.655400 0
M  V30 23 C 7.714400 -8.152900 -4.804200 0
M  V30 24 C 6.568900 -10.221200 -5.128800 0
M  V30 25 C 10.029000 -9.397900 -12.101900 0
M  V30 26 C 1.276500 -8.976000 -5.738400 0
M  V30 27 C -0.227000 -9.045000 -7.763200 0
M  V30 28 C 7.559700 -9.492600 -4.500900 0
M  V30 29 H 5.985500 -10.265600 -9.818600 0
M  V30 30 H 3.752400 -6.907800 -9.426900 0
M  V30 31 H 5.375600 -6.742200 -7.656800 0
M  V30 32 H 4.112200 -7.453000 -6.976100 0
M  V30 33 H 5.135300 -9.443400 -12.144500 0
M  V30 34 H 1.721400 -8.211000 -9.778900 0
M  V30 35 H 7.945500 -11.313500 -10.489000 0
M  V30 36 H 7.512000 -7.409500 -11.481100 0
M  V30 37 H 1.788400 -10.439000 -9.649000 0
M  V30 38 H 3.288600 -10.409500 -9.061600 0
M  V30 39 H 7.006900 -6.634700 -5.952500 0
M  V30 40 H 5.052200 -10.129500 -6.494200 0
M  V30 41 H 10.032600 -11.353800 -11.552000 0
M  V30 42 H 9.625200 -7.437100 -12.478200 0
M  V30 43 H 3.316700 -9.041000 -13.303800 0
M  V30 44 H 4.111700 -7.656900 -13.418500 0
M  V30 45 H 2.528400 -7.691500 -13.644700 0
M  V30 46 H 8.387000 -7.645100 -4.369000 0
M  V30 47 H 6.453800 -11.140000 -4.916700 0
M  V30 48 H 10.885000 -9.406700 -12.513600 0
M  V30 49 H 0.930100 -8.104500 -5.452300 0
M  V30 50 H 0.739600 -9.688700 -5.333600 0
M  V30 51 H 2.208700 -9.066700 -5.452300 0
M  V30 52 H -0.649600 -8.200100 -7.502400 0
M  V30 53 H -0.220600 -9.118900 -8.740000 0
M  V30 54 H -0.732200 -9.793200 -7.381700 0
M  V30 55 H 8.129100 -9.908600 -3.864800 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 2 4
M  V30 4 1 2 5
M  V30 5 1 3 6
M  V30 6 1 3 7
M  V30 7 1 4 8
M  V30 8 1 4 9
M  V30 9 1 5 10
M  V30 10 2 6 11
M  V30 11 1 6 12
M  V30 12 1 8 13
M  V30 13 1 9 14
M  V30 14 1 9 15
M  V30 15 2 10 16
M  V30 16 1 10 17
M  V30 17 1 11 18
M  V30 18 2 12 19
M  V30 19 1 13 20
M  V30 20 1 14 21
M  V30 21 1 15 22
M  V30 22 1 16 23
M  V30 23 2 17 24
M  V30 24 2 18 25
M  V30 25 1 21 26
M  V30 26 1 21 27
M  V30 27 2 23 28
M  V30 28 2 7 8
M  V30 29 1 19 25
M  V30 30 1 21 22
M  V30 31 1 24 28
M  V30 32 1 3 29 CFG=1
M  V30 33 1 4 30 CFG=3
M  V30 34 1 5 31
M  V30 35 1 5 32
M  V30 36 1 7 33
M  V30 37 1 9 34 CFG=3
M  V30 38 1 11 35
M  V30 39 1 12 36
M  V30 40 1 15 37
M  V30 41 1 15 38
M  V30 42 1 16 39
M  V30 43 1 17 40
M  V30 44 1 18 41
M  V30 45 1 19 42
M  V30 46 1 20 43
M  V30 47 1 20 44
M  V30 48 1 20 45
M  V30 49 1 23 46
M  V30 50 1 24 47
M  V30 51 1 25 48
M  V30 52 1 26 49
M  V30 53 1 26 50
M  V30 54 1 26 51
M  V30 55 1 27 52
M  V30 56 1 27 53
M  V30 57 1 27 54
M  V30 58 1 28 55
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(3 3 4 9)
M  V30 END COLLECTION
M  V30 END CTAB
M  END)CTAB";
  v2::FileParsers::MolFileParserParams params;
  params.removeHs = false;
  auto mol = v2::FileParsers::MolFromMolBlock(ctab, params);
  REQUIRE(mol);

  std::vector<DGeomHelpers::UserProvidedTorsion> uts = {
      {
          {0, 1, 4, 9},
          {0.38559172089866883, 0.6641386762221574, 0.8780462436883112,
           0.6808812939306702, 0.7124088530209258, 0.5483229223076913},
          {-1, 1, 1, -1, -1, -1},
      },
      {
          {0, 2, 5, 10},
          {0.5271847026859319, 0.5164502744678869, 0.5676532115509233,
           0.6997359346915658, 0.5365931196993329, 0.6477716042158909},
          {1, -1, 1, 1, -1, 1},
      },
      {
          {0, 2, 5, 11},
          {0.5271847026859319, 0.5164502744678869, 0.5676532115509233,
           0.6997359346915658, 0.5365931196993329, 0.6477716042158909},
          {1, -1, 1, 1, -1, 1},
      },
      {
          {1, 3, 8, 13},
          {0.33571461473723224, 0.6086977355991332, 1.4167321227398417,
           0.6456525541692106, 0.5284760726201473, 0.8024078255039075},
          {-1, 1, 1, -1, -1, -1},
      },
      {
          {1, 4, 9, 15},
          {0.5574253313327392, 0.5976390907251766, 0.5820648760614179,
           0.5859572066810051, 0.530606267231465, 0.5309591327658433},
          {1, 1, 1, 1, -1, 1},
      },
      {
          {1, 4, 9, 16},
          {0.5574253313327392, 0.5976390907251766, 0.5820648760614179,
           0.5859572066810051, 0.530606267231465, 0.5309591327658433},
          {1, 1, 1, 1, -1, 1},
      },
      {
          {3, 7, 12, 19},
          {1.292276352541412, 1.4945356976895072, 1.288246397654784,
           1.3260735515525908, 1.2353931528923832, 1.1559380294779893},
          {1, -1, 1, -1, 1, -1},
      },
  };

  // for (auto &ut : uts) {
  //   for (auto &v : ut.V) {
  //     v *= 0.001;
  //   }
  // }

  DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
  ps.randomSeed = 0xA700F + 1;
  ps.useRandomCoords = false;
  ps.numThreads = 1;
  ps.trackFailures = true;
  ps.verbose = true;
  ps.explicitTorsions = uts;
  ps.maxIterations = 1;
  auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
  std::cerr << "FAILURES: ";
  std::copy(ps.failures.begin(), ps.failures.end(),
            std::ostream_iterator<unsigned int>(std::cerr, " "));
  std::cerr << std::endl;
  REQUIRE(cid == 0);
  std::ofstream("mol.blah4.mol") << MolToV3KMolBlock(*mol) << std::endl;
  // CHECK(cid.size() == 1);
}

TEST_CASE("WOT2") {
  auto mol = "c1ccccc1c1ccccc1"_smiles;
  REQUIRE(mol);
  MolOps::addHs(*mol);

  std::vector<DGeomHelpers::UserProvidedTorsion> uts = {
      {
          {0, 5, 6, 7},
          {0.5574253313327392, 0.5976390907251766, 0.5820648760614179,
           0.5859572066810051, 0.530606267231465, 0.5309591327658433},
          {1, 1, 1, 1, -1, 1},
      },
  };

  // for (auto &ut : uts) {
  //   for (auto &v : ut.V) {
  //     v *= 0.001;
  //   }
  // }

  DGeomHelpers::EmbedParameters ps = DGeomHelpers::ETKDGv3;
  ps.randomSeed = 0xA700F + 1;
  ps.useRandomCoords = false;
  ps.numThreads = 1;
  ps.trackFailures = true;
  ps.verbose = true;
  ps.explicitTorsions = uts;
  ps.maxIterations = 1;
  auto cid = DGeomHelpers::EmbedMolecule(*mol, ps);
  std::cerr << "FAILURES: ";
  std::copy(ps.failures.begin(), ps.failures.end(),
            std::ostream_iterator<unsigned int>(std::cerr, " "));
  std::cerr << std::endl;
  REQUIRE(cid == 0);
  std::ofstream("blah4.mol") << MolToV3KMolBlock(*mol) << std::endl;
  // CHECK(cid.size() == 1);
}
