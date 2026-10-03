//
//  Copyright (C) 2025 RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <catch2/catch_all.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/AvalonGenerator.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <DataStructs/BitOps.h>

using namespace RDKit;

TEST_CASE("Avalon bit and count fingerprints") {
  auto m1 = "c1ccccc1C(=O)NCC(Cl)C"_smiles;
  auto m2 = "CCCCCCCC"_smiles;
  REQUIRE(m1);
  REQUIRE(m2);
  std::unique_ptr<FingerprintGenerator<std::uint32_t>> gen(
      AvalonFP::getAvalonGenerator<std::uint32_t>());

  SECTION("bits") {
    std::unique_ptr<ExplicitBitVect> fp1(gen->getFingerprint(*m1));
    std::unique_ptr<ExplicitBitVect> fp1b(gen->getFingerprint(*m1));
    std::unique_ptr<ExplicitBitVect> fp2(gen->getFingerprint(*m2));
    CHECK(fp1->getNumBits() == 512);
    CHECK(fp1->getNumOnBits() > 20);
    CHECK(*fp1 == *fp1b);
    CHECK(*fp1 != *fp2);
    CHECK(TanimotoSimilarity(*fp1, *fp1) == 1.0);
    CHECK(TanimotoSimilarity(*fp1, *fp2) < 1.0);
  }
  SECTION("counts are consistent with bits") {
    auto fp = gen->getFingerprint(*m1);
    auto cfp = gen->getCountFingerprint(*m1);
    auto counts = AvalonFP::getAvalonCounts(*m1, 512);
    REQUIRE(counts.size() == 512);
    unsigned int nOn = 0;
    for (unsigned int i = 0; i < 512; ++i) {
      CHECK(cfp->getVal(i) == static_cast<int>(counts[i]));
      CHECK(fp->getBit(i) == (counts[i] > 0));
      nOn += counts[i] > 0;
    }
    CHECK(nOn == fp->getNumOnBits());
    bool anyMulti = false;
    for (auto c : counts) {
      anyMulti |= c > 1;
    }
    CHECK(anyMulti);
    delete fp;
    delete cfp;
  }
  SECTION("fingerprint size and flags") {
    std::unique_ptr<FingerprintGenerator<std::uint32_t>> gen2(
        AvalonFP::getAvalonGenerator<std::uint32_t>(1024));
    std::unique_ptr<ExplicitBitVect> fp(gen2->getFingerprint(*m1));
    CHECK(fp->getNumBits() == 1024);
    std::unique_ptr<FingerprintGenerator<std::uint32_t>> gen3(
        AvalonFP::getAvalonGenerator<std::uint32_t>(512, 0x10));  // atom count
    std::unique_ptr<ExplicitBitVect> fp3(gen3->getFingerprint(*m1));
    std::unique_ptr<ExplicitBitVect> fpAll(gen->getFingerprint(*m1));
    CHECK(fp3->getNumOnBits() < fpAll->getNumOnBits());
    CHECK(fp3->getNumOnBits() > 0);
  }
  SECTION("explicit Hs and non-sanitized input") {
    auto mh = "CC(=O)O"_smiles;
    std::unique_ptr<ROMol> mhs(MolOps::addHs(static_cast<const ROMol &>(*mh)));
    std::unique_ptr<ExplicitBitVect> a(gen->getFingerprint(*mh));
    std::unique_ptr<ExplicitBitVect> b(gen->getFingerprint(*mhs));
    CHECK(*a == *b);
    SmilesParserParams ps;
    ps.sanitize = false;
    std::unique_ptr<RWMol> unsan(SmilesToMol("CC(=O)O", ps));
    std::unique_ptr<ExplicitBitVect> c(gen->getFingerprint(*unsan));
    CHECK(*a == *c);
  }
  SECTION("sparse count fingerprint, 64 bit and JSON") {
    std::unique_ptr<FingerprintGenerator<std::uint64_t>> gen64(
        AvalonFP::getAvalonGenerator<std::uint64_t>());
    auto sfp = gen64->getSparseCountFingerprint(*m1);
    auto cnts = AvalonFP::getAvalonCounts(*m1);
    for (unsigned int i = 0; i < cnts.size(); ++i) {
      CHECK(sfp->getVal(i) == static_cast<int>(cnts[i]));
    }
    auto json = generatorToJSON(*gen64);
    auto gen64b = generatorFromJSON(json);
    REQUIRE(gen64b);
    auto sfp2 = gen64b->getSparseCountFingerprint(*m1);
    CHECK(*sfp == *sfp2);
  }
  SECTION("query mode and empty molecules") {
    auto q = AvalonFP::getAvalonCounts(*m1, 512, 0x7FFF, true);
    auto nq = AvalonFP::getAvalonCounts(*m1, 512, 0x7FFF, false);
    CHECK(q != nq);
    ROMol empty;
    std::unique_ptr<ExplicitBitVect> fp(gen->getFingerprint(empty));
    CHECK(fp->getNumOnBits() == 0);
  }
  SECTION("focus atom") {
    auto all = AvalonFP::getAvalonCounts(*m1, 512);
    auto focus = AvalonFP::getAvalonCounts(*m1, 512, 0x7FFF, false, 0);
    CHECK(all != focus);
  }
}
