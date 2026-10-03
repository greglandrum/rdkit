//
//  Copyright (C) 2025 RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Compares the pure-RDKit Avalon fingerprint generator with the original
// implementation from the Avalon toolkit.
#include <catch2/catch_all.hpp>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/AvalonGenerator.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include "AvalonTools.h"

using namespace RDKit;

namespace {
const std::vector<std::string> molSmis = {
    "CCCCCCCC",
    "CC(C)C(=O)O",
    "c1ccccc1C(=O)NCC(Cl)C",
    "CC(=O)Oc1ccccc1C(=O)O",
    "CN1CCC[C@H]1c2cccnc2",
    "C1CC2CCC1CC2",
    "OC(=O)CCc1c[nH]c2ccccc12",
    "[Na+].[O-]C(=O)C",
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "c1ccc2ccccc2c1",
    "c1ccncc1",
    "c1cc[nH]c1",
    "c1ccoc1",
    "c1ccsc1",
    "C1CCCCC1",
    "C1=CC=CC=C1",
    "CC(=O)NC1=CC=C(O)C=C1",
    "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "CS(=O)(=O)Nc1ccccc1",
    "C[N+](C)(C)CC(=O)[O-]",
    "N#CC(C)(C)N=NC(C)(C)C#N",
    "CC1=CC(=O)C=CC1=O",
    "Clc1ccc(cc1)C(c1ccc(Cl)cc1)C(Cl)(Cl)Cl",
    "O=C1CCC(=O)N1",
    "c1ccc2[nH]ccc2c1",
    "c1ccc2c(c1)[nH]c1ccccc12",
    "C1=CC2=CC=CC=C2C=C1",
    "OC(=O)C1=CC=CN=C1",
    "CCOC(=O)C1=C(C)NC(C)=C(C1c1ccccc1[N+]([O-])=O)C(=O)OC",
    "[2H]C([2H])([2H])O",
    "C[C@H](N)C(=O)O",
    "C/C=C/C=O",
    "FC(F)(F)c1ccccc1",
    "C1CC1C(=O)N1CCOCC1",
    "CC(C)(C)c1cc(O)ccc1",
};

std::vector<int> avalonCounts(const ROMol &mol, unsigned int nBits,
                              unsigned int flags) {
  SparseIntVect<std::uint32_t> counts(nBits);
  AvalonTools::getAvalonCountFP(mol, counts, nBits, false, true, flags);
  std::vector<int> res(nBits, 0);
  for (const auto &pr : counts.getNonzeroElements()) {
    res[pr.first] = pr.second;
  }
  return res;
}
}  // namespace

TEST_CASE("Avalon generator bit fingerprints match the Avalon toolkit") {
  for (unsigned int nBits : {512u, 1024u}) {
    for (unsigned int flags :
         {static_cast<unsigned int>(AvalonTools::avalonSSSBits),
          static_cast<unsigned int>(AvalonTools::avalonSimilarityBits)}) {
      std::unique_ptr<FingerprintGenerator<std::uint32_t>> gen(
          AvalonFP::getAvalonGenerator<std::uint32_t>(nBits, flags));
      for (const auto &smi : molSmis) {
        INFO(smi << " nBits=" << nBits << " flags=" << flags);
        std::unique_ptr<ROMol> mol(SmilesToMol(smi));
        REQUIRE(mol);
        ExplicitBitVect ref(nBits);
        AvalonTools::getAvalonFP(*mol, ref, nBits, false, true, flags);
        std::unique_ptr<ExplicitBitVect> fp(gen->getFingerprint(*mol));
        CHECK(*fp == ref);
      }
    }
  }
}

TEST_CASE("Avalon generator count fingerprints match the Avalon toolkit") {
  const unsigned int nBits = 512;
  for (unsigned int flags :
       {static_cast<unsigned int>(AvalonTools::avalonSSSBits),
        static_cast<unsigned int>(AvalonTools::avalonSimilarityBits)}) {
    // the count wrapper in AvalonTools runs its accumulation pass in query mode
    std::unique_ptr<FingerprintGenerator<std::uint32_t>> gen(
        AvalonFP::getAvalonGenerator<std::uint32_t>(nBits, flags, false, true));
    for (const auto &smi : molSmis) {
      INFO(smi << " flags=" << flags);
      std::unique_ptr<ROMol> mol(SmilesToMol(smi));
      REQUIRE(mol);
      auto ref = avalonCounts(*mol, nBits, flags);
      std::unique_ptr<SparseIntVect<std::uint32_t>> fp(
          gen->getCountFingerprint(*mol));
      for (unsigned int i = 0; i < nBits; ++i) {
        INFO("bit " << i);
        CHECK(fp->getVal(i) == ref[i]);
      }
    }
  }
}

TEST_CASE("Avalon generator query fingerprints match the Avalon toolkit") {
  const unsigned int nBits = 512;
  const unsigned int flags = AvalonTools::avalonSSSBits;
  std::unique_ptr<FingerprintGenerator<std::uint32_t>> gen(
      AvalonFP::getAvalonGenerator<std::uint32_t>(nBits, flags, true));
  for (const auto &smi : molSmis) {
    INFO(smi);
    std::unique_ptr<ROMol> mol(SmilesToMol(smi));
    REQUIRE(mol);
    ExplicitBitVect ref(nBits);
    AvalonTools::getAvalonFP(*mol, ref, nBits, true, true, flags);
    std::unique_ptr<ExplicitBitVect> fp(gen->getFingerprint(*mol));
    CHECK(*fp == ref);

    SparseIntVect<std::uint32_t> refCounts(nBits);
    AvalonTools::getAvalonCountFP(*mol, refCounts, nBits, true, true, flags);
    std::unique_ptr<SparseIntVect<std::uint32_t>> cfp(
        gen->getCountFingerprint(*mol));
    for (unsigned int i = 0; i < nBits; ++i) {
      CHECK(cfp->getVal(i) == refCounts.getVal(i));
    }
  }
}
