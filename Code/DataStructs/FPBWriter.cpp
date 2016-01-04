//
// Copyright (c) 2016 Greg Landrum
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Implementation details here are taken from the file fpb_io.py from chemfp
// (www.chemfp.org)
// Many thanks to Andrew Dalke for creating such great software and for
// helping explain the FPB implementation

#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>

#include <RDGeneral/Invariant.h>
#include <RDGeneral/StreamOps.h>
#include <RDGeneral/Ranking.h>
#include "FPBWriter.h"

namespace RDKit {

namespace detail {
const unsigned int magicSize = 8;
const std::string FPB_MAGIC("FPB1\r\n\0\0", 8);
const unsigned int tagNameSize = 4;

struct FPBWriter_impl {
  unsigned int len;
  unsigned int nBits;
  boost::uint32_t numBytesStoredPerFingerprint;
  std::vector<boost::uint32_t> popCountOffsets;
  boost::uint32_t num4ByteElements, num8ByteElements;  // for finding ids
};

// caller is responsible for delete []'ing the result
boost::uint8_t *bitsetToBytes(const boost::dynamic_bitset<> &bitset) {
  unsigned int nBits = bitset.size();
  unsigned int nBytes = nBits / 8;
  unsigned int nBlocks = nBytes / sizeof(boost::dynamic_bitset<>::block_type);

  boost::uint8_t *res = new boost::uint8_t[nBytes];
  boost::to_block_range(bitset, (boost::dynamic_bitset<>::block_type *)res);
  return res;
}

}  // end of detail namespace

void FPBWriter::init() {
  PRECONDITION(dp_strm, "no stream");
  dp_impl = new detail::FPBWriter_impl;

  char magic[detail::magicSize];
  df_init = true;
};

void FPBWriter::destroy() {
  if (dp_impl) {
  }
  delete dp_impl;
};

}  // end of RDKit namespace
