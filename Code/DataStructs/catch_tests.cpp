//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <catch2/catch_all.hpp>

#include <RDGeneral/test.h>
#include "BitVects.h"
#include "BitOps.h"
#include "BitVectUtils.h"
#include "SparseIntVect.h"
#include <limits>

using namespace RDKit;

TEST_CASE("special cases for the limits of sparse vectors") {
  SECTION("SparseBitVect") {
    SparseBitVect sbv(std::numeric_limits<unsigned int>::max());
    CHECK(sbv.getNumBits() == std::numeric_limits<unsigned int>::max());
    CHECK(!sbv.setBit(std::numeric_limits<unsigned int>::max()));
    CHECK(sbv.getBit(std::numeric_limits<unsigned int>::max()) == 1);
  }
}

TEST_CASE("getOnBits") {
  SECTION("SparseBitVect") {
    SparseBitVect bv(100);
    bv.setBit(0);
    bv.setBit(12);
    bv.setBit(17);

    IntVect obl1;
    bv.getOnBits(obl1);
    CHECK(obl1.size() == 3);
    auto obl2 = bv.getOnBits();
    CHECK(obl2.size() == 3);
    CHECK(obl1 == obl2);
  }
  SECTION("ExplicitBitVect") {
    ExplicitBitVect bv(100);
    bv.setBit(0);
    bv.setBit(12);
    bv.setBit(17);

    IntVect obl1;
    bv.getOnBits(obl1);
    CHECK(obl1.size() == 3);
    auto obl2 = bv.getOnBits();
    CHECK(obl2.size() == 3);
    CHECK(obl1 == obl2);
  }
}