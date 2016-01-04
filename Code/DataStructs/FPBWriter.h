//
// Copyright (c) 2016 Greg Landrum
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_FPBWRITER_H_JAN2016
#define RD_FPBWRITER_H_JAN2016

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <RDGeneral/BadFileException.h>
#include <DataStructs/ExplicitBitVect.h>

#include <boost/cstdint.hpp>

namespace RDKit {
namespace detail {
struct FPBWriter_impl;
}
class FPBWriter {
 public:
  FPBWriter() : dp_strm(NULL), dp_impl(NULL), df_owner(false), df_init(false){};
  FPBWriter(const char *fname) { _initFromFilename(fname); };
  FPBWriter(const std::string &fname) { _initFromFilename(fname.c_str()); };
  FPBWriter(std::ostream *outStream, bool takeOwnership = true)
      : dp_strm(outStream), df_owner(takeOwnership), df_init(false){};
  ~FPBWriter() {
    destroy();
    if (df_owner) delete dp_strm;
    dp_strm = NULL;
    df_init = false;
  };
  void init();

 private:
  std::ostream *dp_strm;
  detail::FPBWriter_impl *dp_impl;  // implementation details
  bool df_owner;
  bool df_init;

  // disable automatic copy constructors and assignment operators
  // for this class and its subclasses.  They will likely be
  // carrying around stream pointers and copying those is a recipe
  // for disaster.
  FPBWriter(const FPBWriter &);
  FPBWriter &operator=(const FPBWriter &);
  void destroy();
  void _initFromFilename(const char *fname) {
    std::ostream *tmpStream = static_cast<std::ostream *>(new std::ofstream(
        fname,
        std::ios_base::binary | std::ios_base::out | std::ios_base::trunc));
    if (!tmpStream || (!(*tmpStream)) || (tmpStream->bad())) {
      std::ostringstream errout;
      errout << "Bad input file " << fname;
      throw BadFileException(errout.str());
    }
    dp_strm = tmpStream;
    df_owner = true;
    df_init = false;
  }
};
}
#endif
