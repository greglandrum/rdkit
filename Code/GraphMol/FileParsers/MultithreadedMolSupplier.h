//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifdef RDK_THREADSAFE_SSS
#ifndef MULTITHREADED_MOL_SUPPLIER
#define MULTITHREADED_MOL_SUPPLIER

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/ConcurrentQueue.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/RDThreads.h>
#include <RDGeneral/StreamOps.h>

#include <boost/tokenizer.hpp>
#include <cstdlib>
#include <functional>
#include <memory>
#include <sstream>

#include "FileParsers.h"
#include "MolSupplier.h"

typedef boost::tokenizer<boost::char_separator<char>> tokenizer;

namespace RDKit {
class RDKIT_FILEPARSERS_EXPORT MultithreadedMolSupplier : public MolSupplier {
  // this is an abstract base class to concurrently supply molecules one at a
  // time
 public:
  MultithreadedMolSupplier(){};
  virtual ~MultithreadedMolSupplier(){};
  //! reads lines from input stream to populate the input queue
  void inputProducer();
  //! parses lines from the input queue converting them to ROMol objects
  //! populating the output queue
  void inputConsumer(size_t id);
  //! pop elements from the output queue
  ROMol *next();
  //! starts reader and writer threads
  void startThreads();
  //! use atEnd method of the
  bool atEnd();

 private:
  // disable automatic copy constructors and assignment operators
  // for this class and its subclasses.  They will likely be
  // carrying around stream pointers and copying those is a recipe
  // for disaster.
  MultithreadedMolSupplier(const MultithreadedMolSupplier &);
  MultithreadedMolSupplier &operator=(const MultithreadedMolSupplier &);

 private:
  virtual void reset();
  virtual void init() = 0;
  virtual bool getEnd() = 0;
  virtual bool extractNextRecord(std::string &record,
                                 unsigned int &lineNum) = 0;
  virtual ROMol *processMoleculeRecord(const std::string &record,
                                       unsigned int lineNum) = 0;

 protected:
  const int d_numReaderThread = 1;  // fix number of reader threads to 1
  int d_numWriterThreads;           // number of writer threads
  size_t d_sizeInputQueue;          // size of input queue
  size_t d_sizeOutputQueue;         // size of output queue
  ConcurrentQueue<std::tuple<std::string, unsigned int>>
      *d_inputQueue;  // concurrent input queue
  ConcurrentQueue<std::shared_ptr<ROMol>>
      *d_outputQueue;  // concurrent output queue
};
}  // namespace RDKit
#endif
#endif
