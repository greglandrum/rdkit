#pragma once

#include <GraphMol/FileParsers/MolSupplier.h>

namespace RDKit {
class MAESupplier : public MolSupplier {
  // this is an abstract base class to supply molecules one at a time
 public:
  MAESupplier(){};
  virtual ~MAESupplier(){};
  virtual void init() = 0;
  virtual void reset() = 0;
  virtual bool atEnd() = 0;
  virtual ROMol *next() = 0;

 private:
 protected:
};

}  // end of namespace
