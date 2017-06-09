if(RDK_BUILD_MAESTRO_SUPPORT)
if(NOT DEFINED MAEPARSER_DIR)
  set(MAEPARSER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/maeparser")
  set(fileToCheck "${MAEPARSER_DIR}/MaeParser.hpp")
  set(needDownload "TRUE")
  if(EXISTS "${fileToCheck}")
      set(needDownload "FALSE")
  endif()
else()
  string(REGEX REPLACE "\\\\" "/" MAEPARSER_DIR ${MAEPARSER_DIR})
  set(needDownload "FALSE")
endif()

if(needDownload)
  # FIX: This is completely untested and likely does not work
  if(NOT DEFINED MAEPARSER_URL)
    set(MAEPARSER_URL "https://codeload.github.com/schrodinger/maeparser/zip/master")
  endif()
  if(NOT DEFINED MAEPARSER_MD5SUM)
    set(MAEPARSER_MD5SUM "3f6f16010a7e8d1b94516ff20886d965")
  endif()
  downloadAndCheckMD5(${MAEPARSER_URL} ${MAEPARSER_DIR} ${MAEPARSER_MD5SUM})
  execute_process(COMMAND ${CMAKE_COMMAND} -E unzip
    ${CMAKE_CURRENT_SOURCE_DIR}/maeparser-master.zip
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

add_subdirectory(maeparser)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${MAEPARSER_DIR})

rdkit_headers(MAESupplier.h DEST GraphMol)
#rdkit_test(testAvalonLib1 test1.cpp
#           LINK_LIBRARIES AvalonLib avalon_clib FileParsers SmilesParse GraphMol DataStructs RDGeometryLib RDGeneral)


endif(RDK_BUILD_MAESTRO_SUPPORT)
