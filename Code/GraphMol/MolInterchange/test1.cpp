//
//  Copyright (C) 2018 Greg Landrum
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolInterchange/MolInterchange.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/LocaleSwitcher.h>
#include <clocale>
#include <cstdlib>

#include <string>
#include <fstream>

using namespace RDKit;

void test1() {
  BOOST_LOG(rdInfoLog) << "test1: basics" << std::endl;

  std::string rdbase = getenv("RDBASE");
  std::string fName =
      rdbase + "/Code/GraphMol/MolInterchange/test_data/test1.json";
  std::ifstream inStream(fName);
  if (!inStream || (inStream.bad())) {
    std::ostringstream errout;
    errout << "Bad input file " << fName;
    throw BadFileException(errout.str());
  }

  std::vector<boost::shared_ptr<RWMol>> mols =
      MolInterchange::JSONDataStreamToMols(&inStream);
  TEST_ASSERT(mols.size() == 1);
  RWMol *m = mols[0].get();
  TEST_ASSERT(m);
  m->debugMol(std::cerr);
  TEST_ASSERT(m->getNumAtoms() == 15);

  BOOST_LOG(rdInfoLog) << "done" << std::endl;
}

void RunTests() {
#if 1
  test1();
#endif
}

#if 0
// must be in German Locale for test...
void testLocaleSwitcher() {
  float d = -1.0;
  char buffer[1024];
  sprintf(buffer, "%0.2f", d);
  if (std::string(buffer) != "-1,00") {
    BOOST_LOG(rdInfoLog) << " ---- no German locale support (skipping) ---- "
                         << std::endl;
    return;
  }

  {
    RDKit::Utils::LocaleSwitcher ls;
    sprintf(buffer, "%0.2f", d);
    CHECK_INVARIANT(std::string(buffer) == "-1.00", "Locale Switcher Fail");
    // test locale switcher recursion
    {
      RDKit::Utils::LocaleSwitcher ls;
      sprintf(buffer, "%0.2f", d);
      CHECK_INVARIANT(std::string(buffer) == "-1.00", "Locale Switcher Fail");
    }
    // should still be in the "C" variant
    sprintf(buffer, "%0.2f", d);
    CHECK_INVARIANT(std::string(buffer) == "-1.00", "Locale Switcher Fail");
  }

  // Should be back in German Locale
  sprintf(buffer, "%0.2f", d);
  CHECK_INVARIANT(std::string(buffer) == "-1,00", "Locale Switcher Fail");
}

#ifdef RDK_TEST_MULTITHREADED

#include <RDGeneral/BoostStartInclude.h>
#include <boost/thread.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace {
void runblock() { testLocaleSwitcher(); }
}
void testMultiThreadedSwitcher() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test multithreading Locale Switching"
                        << std::endl;

  boost::thread_group tg;
  unsigned int count = 100;
  for (unsigned int i = 0; i < count; ++i) {
    tg.add_thread(new boost::thread(runblock));
  }
  tg.join_all();
  BOOST_LOG(rdErrorLog) << "    Test multithreading (Done)" << std::endl;
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
}

#else

void testMultiThreadedSwitcher() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdInfoLog) << " ---- Multithreaded tests disabled ---- "
                       << std::endl;
}
#endif
#endif

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  RDLog::InitLogs();
  RunTests();  // run with C locale

  return 0;
}
