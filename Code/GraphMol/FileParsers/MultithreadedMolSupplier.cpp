//
//  Copyright (C) 2020 Shrey Aryan
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MultithreadedMolSupplier.h"
namespace RDKit {
void MultithreadedMolSupplier::reader() {
  std::string record;
  unsigned int lineNum, index;
  while (extractNextRecord(record, lineNum, index)) {
    auto r = std::tuple<std::string, unsigned int, unsigned int>{
        record, lineNum, index};
    d_inputQueue->push(r);
  }
  d_inputQueue->setDone();
}

void MultithreadedMolSupplier::writer() {
  std::tuple<std::string, unsigned int, unsigned int> r;
  while (d_inputQueue->pop(r)) {
    ROMol* mol = processMoleculeRecord(std::get<0>(r), std::get<1>(r));
    auto temp = std::tuple<ROMol*, std::string, unsigned int>{
        mol, std::get<0>(r), std::get<2>(r)};
    d_outputQueue->push(temp);
  }

  if (d_threadCounter != d_numWriterThreads) {
    ++d_threadCounter;
  } else {
    d_outputQueue->setDone();
  }
}

ROMol* MultithreadedMolSupplier::next() {
  std::tuple<ROMol*, std::string, unsigned int> r;
  if (d_outputQueue->pop(r)) {
    ROMol* mol = std::get<0>(r);
    d_lastItemText = std::get<1>(r);
    d_lastRecordId = std::get<2>(r);
    return mol;
  }
  return nullptr;
}

void MultithreadedMolSupplier::endThreads() {
  d_readerThread.join();
  for (auto& thread : d_writerThreads) {
    thread.join();
  }
}

void MultithreadedMolSupplier::startThreads() {
  //! run the reader function in a seperate thread
  d_readerThread = std::thread(&MultithreadedMolSupplier::reader, this);
  //! run the writer function in seperate threads
  for (unsigned int i = 0; i < d_numWriterThreads; i++) {
    d_writerThreads.emplace_back(
        std::thread(&MultithreadedMolSupplier::writer, this));
  }
}

bool MultithreadedMolSupplier::atEnd() {
  return (d_outputQueue->isEmpty() && d_outputQueue->getDone());
}

void MultithreadedMolSupplier::reset() { ; }

}  // namespace RDKit