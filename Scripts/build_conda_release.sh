#!/bin/bash
condapath=`which conda`
BINDIR=`dirname $condapath`
PYROOT=`dirname $BINDIR`
if [[ -d builds ]]; then
  echo "reusing existing builds directory"
else
  mkdir builds
fi
cd builds
echo `dirname "${BASH_SOURCE[0]}"` 
exit
SRCDIR=../../
BLDDIR=build_conda_release

if [[ -d $BLDDIR ]]; then
  echo "reusing existing conda build directory"
else
  mkdir $BLDDIR
fi

#echo "Cleaning old directory"
#rm -rf $BLDDIR
#mkdir $BLDDIR
cd $BLDDIR
echo "cmake creating"
INSTALL=`pwd`/rdkit_build
BOOST_ROOT=$PYROOT

if [ "$(uname -s)" = "Darwin" ]; then
  CXX="/usr/bin/c++ -stdlib=libstdc++"
  PYLIB=$PYROOT/lib/libpython3.5m.dylib
else
  PYLIB=$PYROOT/lib/libpython3.5m.dylib
fi

cmake -DPYTHON_EXECUTABLE=$PYROOT/bin/python  \
          -DAVALONTOOLS_DIR=$SRCDIR/External/AvalonTools/distrib/SourceDistribution -DRDK_BUILD_AVALON_SUPPORT=ON \
          -DRDK_BUILD_INCHI_SUPPORT=ON \
          -DRDK_BUILD_CAIRO_SUPPORT=ON -DRDK_BUILD_THREADSAFE_SSS=ON -DRDK_TEST_MULTITHREADED=ON \
	  -DPYTHON_LIBRARY=$PYLIB \
	  -DPYTHON_INCLUDE_DIR=$PYROOT/include/python3.5m \
	  -DPYTHON_NUMPY_INCLUDE_PATH=$PYROOT/lib/python3.5/site-packages/numpy/core/include/ \
	  -DCMAKE_BUILD_TYPE=Release \
	  -DRDK_INSTALL_INTREE=OFF \
	  -DCMAKE_INSTALL_PREFIX=$INSTALL \
          -DBOOST_ROOT=$PYROOT \
          -DBoost_NO_BOOST_CMAKE=FALSE \
          -DBoost_NO_SYSTEM_PATHS=on \
          -DBoost_USE_STATIC_LIBS=off \
	  $SRCDIR && \
    make -j 2 install && \
    RDBASE=`pwd`/$SRCDIR DYLD_FALLBACK_LIBRARY_PATH=`pwd`/rdkit_build/lib:$PYROOT/lib PYTHONPATH=`pwd`/rdkit_build/lib/python3.5/site-packages ctest -j2 \
    echo "To use this environment, do:" && \
    echo "export RDBASE="`pwd`/$SRCDIR && \
    echo "export DYLD_FALLBACK_LIBRARY_PATH="`pwd`/rdkit_build/lib:$PYROOT/lib && \
    echo "export PYTHONPATH=`pwd`/rdkit_build/lib/python3.5/site-packages"
