rdkit_library(DataStructs 
              BitVect.cpp SparseBitVect.cpp ExplicitBitVect.cpp Utils.cpp
              base64.cpp BitOps.cpp DiscreteDistMat.cpp
              DiscreteValueVect.cpp FPBReader.cpp
              LINK_LIBRARIES RDGeneral ${RDKit_THREAD_LIBS})

rdkit_headers(base64.h
              BitOps.h
              BitVect.h
              BitVects.h
              BitVectUtils.h
              DatastructsException.h
              DiscreteDistMat.h
              DiscreteValueVect.h
              ExplicitBitVect.h
              SparseBitVect.h
              SparseIntVect.h
              FPBReader.h DEST DataStructs)

rdkit_test(testDataStructs testDatastructs.cpp 
           LINK_LIBRARIES DataStructs RDGeneral)

rdkit_test(testFPB testFPB.cpp 
           LINK_LIBRARIES DataStructs RDGeneral ${RDKit_THREAD_LIBS})

add_subdirectory(Wrap)
