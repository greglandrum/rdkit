//
//   Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _GRID3D_H_20050124_1113
#define _GRID3D_H_20050124_1113
#include <exception>
#include <string>
#include <utility>

namespace RDKit {
class DiscreteValueVect;
}
namespace RDGeom {
class Point3D;

class RDKIT_RDGEOMETRYLIB_EXPORT GridException : public std::exception {
 public:
  //! construct with an error message
  GridException(const char *msg) : _msg(msg) {}
  //! construct with an error message
  GridException(std::string msg) : _msg(std::move(msg)) {}
  //! get the error message
  const char *what() const noexcept override { return _msg.c_str(); }
  ~GridException() noexcept override = default;

 private:
  std::string _msg;
};

//! Virtual base class for a grid object
template <typename VectorType, typename ValueType1, typename ValueType2>
class RDKIT_RDGEOMETRYLIB_EXPORT Grid3D {
 public:
  virtual ~Grid3D(){};
  virtual int getGridPointIndex(const Point3D &point) const = 0;
  virtual ValueType1 getVal(const Point3D &point) const = 0;
  virtual void setVal(const Point3D &point, ValueType2 val) = 0;

  virtual Point3D getGridPointLoc(unsigned int pointId) const = 0;
  virtual ValueType2 getVal(unsigned int pointId) const = 0;
  virtual void setVal(unsigned int pointId, ValueType2 val) = 0;
  virtual unsigned int getSize() const = 0;
  virtual const VectorType *getOccupancyVect() const = 0;
};
}  // namespace RDGeom

#endif
