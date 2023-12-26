//
//  Copyright (c) 2010-2021 Novartis Institutes for BioMedical Research Inc.
//    and other RDKit contributors
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

// PostgreSQL 14 on Windows uses a hack to redefine the stat struct
// The hack assumes that sys/stat.h will be imported for the first
// time by win32_port.h, which is not necessarily the case
// So we need to set the stage for the hack or it will fail
#ifdef _WIN32
#define fstat microsoft_native_fstat
#define stat microsoft_native_stat
#include <sys/stat.h>
#ifdef __MINGW32__
#ifndef HAVE_GETTIMEOFDAY
#define HAVE_GETTIMEOFDAY 1
#endif
#endif
#endif

#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include <GraphMol/GenericGroups/GenericGroups.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/MACCS.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/MolHash/MolHash.h>
#include <GraphMol/FMCS/FMCS.h>
#include <DataStructs/BitOps.h>
#include <DataStructs/SparseIntVect.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/MolEnumerator/MolEnumerator.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/integer_traits.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <GraphMol/ChemReactions/ReactionFingerprints.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>


// see above comment on the PostgreSQL hack
#ifdef _WIN32
#undef fstat
#undef stat
#endif

#include "rdkit.h"
#include "guc.h"
#include "bitstring.h"

using namespace std;
using namespace RDKit;

constexpr unsigned int pickleForQuery =
    PicklerOps::PropertyPickleOptions::MolProps |
    PicklerOps::PropertyPickleOptions::AtomProps |
    PicklerOps::PropertyPickleOptions::BondProps |
    PicklerOps::PropertyPickleOptions::PrivateProps |
    PicklerOps::PropertyPickleOptions::QueryAtomData;
constexpr unsigned int pickleDefault =
    PicklerOps::PropertyPickleOptions::MolProps |
    PicklerOps::PropertyPickleOptions::PrivateProps;

class ByteA2 : public std::string {
 public:
  ByteA2() : string(){};
  ByteA2(bytea *b) : string(VARDATA(b), VARSIZE(b) - VARHDRSZ){};
  ByteA2(string &s) : string(s){};

  /*
   * Convert string to bytea. Convertaion is in pgsql's memory
   */
  bytea *toByteA2() {
    bytea *res;
    int len;

    len = this->size();
    res = (bytea *)palloc(VARHDRSZ + len);
    memcpy(VARDATA(res), this->data(), len);
    SET_VARSIZE(res, VARHDRSZ + len);

    return res;
  };

  /* Just the copy of string's method */
  ByteA2 &operator=(const string &__str) {
    return (ByteA2 &)this->assign(__str);
  };
};

/*
 * Constant io
 */
static string StringData;

/*
 * Real sparse vector
 */

typedef SparseIntVect<std::uint32_t> SparseFP;



extern "C" char *ReactionGetSVG(CChemicalReaction i, unsigned int w,
                                unsigned int h, bool highlightByReactant,
                                const char *params) {
  auto *rxn = (ChemicalReaction *)i;

  MolDraw2DSVG drawer(w, h);
  if (params && strlen(params)) {
    try {
      MolDraw2DUtils::updateDrawerParamsFromJSON(drawer, params);
    } catch (...) {
      elog(WARNING,
           "adjustQueryProperties: Invalid argument \'params\' ignored");
    }
  }
  drawer.drawReaction(*rxn, highlightByReactant);
  drawer.finishDrawing();
  std::string txt = drawer.getDrawingText();
  return strdup(txt.c_str());
}

/* chemical reactions */

extern "C" void freeChemReaction(CChemicalReaction data) {
  auto *rxn = (ChemicalReaction *)data;
  delete rxn;
}

extern "C" CChemicalReaction constructChemReact(Reaction *data) {
  auto *rxn = new ChemicalReaction();

  try {
    ByteA2 b(data);
    ReactionPickler::reactionFromPickle(b, rxn);
  } catch (ReactionPicklerException &e) {
    elog(ERROR, "reactionFromPickle: %s", e.what());
  } catch (...) {
    elog(ERROR, "constructChemReact: Unknown exception");
  }

  return (CChemicalReaction)rxn;
}

extern "C" Reaction *deconstructChemReact(CChemicalReaction data) {
  auto *rxn = (ChemicalReaction *)data;
  ByteA2 b;

  try {
    ReactionPickler::pickleReaction(rxn, b);
  } catch (ReactionPicklerException &e) {
    elog(ERROR, "pickleReaction: %s", e.what());
  } catch (...) {
    elog(ERROR, "deconstructChemReact: Unknown exception");
  }

  return (Reaction *)b.toByteA2();
}

extern "C" CChemicalReaction parseChemReactText(char *data, bool asSmarts,
                                                bool warnOnFail) {
  ChemicalReaction *rxn = nullptr;

  try {
    if (asSmarts) {
      rxn = RxnSmartsToChemicalReaction(data);
    } else {
      rxn = RxnSmartsToChemicalReaction(data, nullptr, true);
    }
    if (getInitReaction()) {
      rxn->initReactantMatchers();
    }
    if (getMoveUnmappedReactantsToAgents() && hasReactionAtomMapping(*rxn)) {
      rxn->removeUnmappedReactantTemplates(getThresholdUnmappedReactantAtoms());
    }
  } catch (...) {
    rxn = nullptr;
  }
  if (rxn == nullptr) {
    if (warnOnFail) {
      ereport(WARNING,
              (errcode(ERRCODE_WARNING),
               errmsg("could not create chemical reaction from SMILES '%s'",
                      data)));
    } else {
      ereport(ERROR,
              (errcode(ERRCODE_DATA_EXCEPTION),
               errmsg("could not create chemical reaction  from SMILES '%s'",
                      data)));
    }
  }

  return (CChemicalReaction)rxn;
}

extern "C" CChemicalReaction parseChemReactBlob(char *data, int len) {
  ChemicalReaction *rxn = nullptr;

  try {
    string binStr(data, len);
    rxn = new ChemicalReaction(binStr);
    if (getInitReaction()) {
      rxn->initReactantMatchers();
    }
    if (getMoveUnmappedReactantsToAgents() && hasReactionAtomMapping(*rxn)) {
      rxn->removeUnmappedReactantTemplates(getThresholdUnmappedReactantAtoms());
    }
  } catch (...) {
    ereport(ERROR,
            (errcode(ERRCODE_DATA_EXCEPTION),
             errmsg("problem generating chemical reaction from blob data")));
  }
  if (rxn == nullptr) {
    ereport(ERROR, (errcode(ERRCODE_DATA_EXCEPTION),
                    errmsg("blob data could not be parsed")));
  }

  return (CChemicalReaction)rxn;
}

extern "C" char *makeChemReactText(CChemicalReaction data, int *len,
                                   bool asSmarts) {
  auto *rxn = (ChemicalReaction *)data;

  try {
    if (!asSmarts) {
      StringData = ChemicalReactionToRxnSmiles(*rxn);
    } else {
      StringData = ChemicalReactionToRxnSmarts(*rxn);
    }
  } catch (...) {
    ereport(WARNING, (errcode(ERRCODE_WARNING),
                      errmsg("makeChemReactText: problems converting chemical "
                             "reaction  to SMILES/SMARTS")));
    StringData = "";
  }

  *len = StringData.size();
  return (char *)StringData.c_str();
}

extern "C" char *makeChemReactBlob(CChemicalReaction data, int *len) {
  auto *rxn = (ChemicalReaction *)data;
  StringData.clear();
  try {
    ReactionPickler::pickleReaction(*rxn, StringData);
  } catch (...) {
    elog(ERROR, "makeChemReactBlob: Unknown exception");
  }

  *len = StringData.size();
  return (char *)StringData.data();
}

extern "C" CChemicalReaction parseChemReactCTAB(char *data, bool warnOnFail) {
  ChemicalReaction *rxn = nullptr;

  try {
    rxn = RxnBlockToChemicalReaction(data);
    if (getInitReaction()) {
      rxn->initReactantMatchers();
    }
    if (getMoveUnmappedReactantsToAgents() && hasReactionAtomMapping(*rxn)) {
      rxn->removeUnmappedReactantTemplates(getThresholdUnmappedReactantAtoms());
    }
  } catch (...) {
    rxn = nullptr;
  }
  if (rxn == nullptr) {
    if (warnOnFail) {
      ereport(WARNING,
              (errcode(ERRCODE_WARNING),
               errmsg("could not create reaction from CTAB '%s'", data)));

    } else {
      ereport(ERROR,
              (errcode(ERRCODE_DATA_EXCEPTION),
               errmsg("could not create reaction from CTAB '%s'", data)));
    }
  }

  return (CChemicalReaction)rxn;
}

extern "C" char *makeCTABChemReact(CChemicalReaction data, int *len) {
  auto *rxn = (ChemicalReaction *)data;

  try {
    StringData = ChemicalReactionToRxnBlock(*rxn);
  } catch (...) {
    ereport(
        WARNING,
        (errcode(ERRCODE_WARNING),
         errmsg("makeCTABChemReact: problems converting reaction to CTAB")));
    StringData = "";
  }

  *len = StringData.size();
  return (char *)StringData.c_str();
}

extern "C" int ChemReactNumReactants(CChemicalReaction crxn) {
  const ChemicalReaction *rxn = (ChemicalReaction *)crxn;
  return rxn->getNumReactantTemplates();
}

extern "C" int ChemReactNumProducts(CChemicalReaction crxn) {
  const ChemicalReaction *rxn = (ChemicalReaction *)crxn;
  return rxn->getNumProductTemplates();
}

extern "C" int ChemReactNumAgents(CChemicalReaction crxn) {
  const ChemicalReaction *rxn = (ChemicalReaction *)crxn;
  return rxn->getNumAgentTemplates();
}

extern "C" bytea *makeReactionSign(CChemicalReaction data) {
  auto *rxn = (ChemicalReaction *)data;
  ExplicitBitVect *res = nullptr;
  bytea *ret = nullptr;

  try {
    RDKit::ReactionFingerprintParams params;
    params.fpType = static_cast<FingerprintType>(getReactionSubstructFpType());
    params.fpSize = getReactionSubstructFpSize();
    params.includeAgents = (!getIgnoreReactionAgents());
    params.bitRatioAgents = getReactionStructuralFPAgentBitRatio();
    res = RDKit::StructuralFingerprintChemReaction(*rxn, params);

    if (res) {
      std::string sres = BitVectToBinaryText(*res);

      unsigned int varsize = VARHDRSZ + sres.size();
      ret = (bytea *)palloc0(varsize);
      memcpy(VARDATA(ret), sres.data(), sres.size());
      SET_VARSIZE(ret, varsize);

      delete res;
      res = nullptr;
    }
  } catch (...) {
    elog(ERROR, "makeReactionSign: Unknown exception");
    if (res) {
      delete res;
    }
  }
  return ret;
}

extern "C" int ReactionSubstruct(CChemicalReaction rxn,
                                 CChemicalReaction rxn2) {
  auto *rxnm = (ChemicalReaction *)rxn;
  auto *rxn2m = (ChemicalReaction *)rxn2;

  /* Reaction search */
  if (rxn2m->getNumReactantTemplates() != 0 &&
      rxn2m->getNumProductTemplates() != 0) {
    return hasReactionSubstructMatch(*rxnm, *rxn2m,
                                     (!getIgnoreReactionAgents()));
  }
  /* Product search */
  if (rxn2m->getNumReactantTemplates() == 0 &&
      rxn2m->getNumProductTemplates() != 0) {
    if (rxn2m->getNumAgentTemplates() != 0 && !getIgnoreReactionAgents()) {
      return (hasProductTemplateSubstructMatch(*rxnm, *rxn2m) &&
              hasAgentTemplateSubstructMatch(*rxnm, *rxn2m));
    }
    return hasProductTemplateSubstructMatch(*rxnm, *rxn2m);
  }
  /* Reactant search */
  if (rxn2m->getNumReactantTemplates() != 0 &&
      rxn2m->getNumProductTemplates() == 0) {
    if (rxn2m->getNumAgentTemplates() != 0 && !getIgnoreReactionAgents()) {
      return (hasReactantTemplateSubstructMatch(*rxnm, *rxn2m) &&
              hasAgentTemplateSubstructMatch(*rxnm, *rxn2m));
    }
    return hasReactantTemplateSubstructMatch(*rxnm, *rxn2m);
  }
  /* Agent search */
  if (rxn2m->getNumReactantTemplates() == 0 &&
      rxn2m->getNumProductTemplates() == 0 &&
      rxn2m->getNumAgentTemplates() != 0) {
    return hasAgentTemplateSubstructMatch(*rxnm, *rxn2m);
  }

  return false;
}

extern "C" int ReactionSubstructFP(CChemicalReaction rxn,
                                   CChemicalReaction rxnquery) {
  auto *rxnm = (ChemicalReaction *)rxn;
  auto *rxnqm = (ChemicalReaction *)rxnquery;

  RDKit::ReactionFingerprintParams params;
  params.fpType = static_cast<FingerprintType>(getReactionSubstructFpType());
  params.fpSize = getReactionSubstructFpSize();
  params.includeAgents = (!getIgnoreReactionAgents());
  params.bitRatioAgents = getReactionStructuralFPAgentBitRatio();

  ExplicitBitVect *fp1 = StructuralFingerprintChemReaction(*rxnm, params);
  ExplicitBitVect *fp2 = StructuralFingerprintChemReaction(*rxnqm, params);

  if (fp1->getNumOnBits() < fp2->getNumOnBits()) {
    return false;
  }
  for (unsigned i = 0; i < fp1->getNumBits(); i++) {
    if ((fp1->getBit(i) & fp2->getBit(i)) != fp2->getBit(i)) {
      return false;
    }
  }
  return true;
}

// some helper functions in anonymous namespace
namespace {

struct MoleculeDescriptors {
  MoleculeDescriptors() {}
  unsigned nAtoms{0};
  unsigned nBonds{0};
  unsigned nRings{0};
  double MW{0.0};
};

MoleculeDescriptors *calcMolecularDescriptorsReaction(
    RDKit::ChemicalReaction *rxn, RDKit::ReactionMoleculeType t) {
  auto *des = new MoleculeDescriptors();
  auto begin = getStartIterator(*rxn, t);
  auto end = getEndIterator(*rxn, t);
  for (; begin != end; ++begin) {
    des->nAtoms += begin->get()->getNumHeavyAtoms();
    des->nBonds += begin->get()->getNumBonds(true);
    des->MW = RDKit::Descriptors::calcAMW(*begin->get(), true);
    if (!begin->get()->getRingInfo()->isSssrOrBetter()) {
      begin->get()->updatePropertyCache();
      RDKit::MolOps::findSSSR(*begin->get());
    }
    des->nRings += begin->get()->getRingInfo()->numRings();
  }
  return des;
}

int compareMolDescriptors(const MoleculeDescriptors &md1,
                          const MoleculeDescriptors &md2) {
  int res = md1.nAtoms - md2.nAtoms;
  if (res) {
    return res;
  }
  res = md1.nBonds - md2.nBonds;
  if (res) {
    return res;
  }
  res = md1.nRings - md2.nRings;
  if (res) {
    return res;
  }
  res = int(md1.MW - md2.MW);
  if (res) {
    return res;
  }
  return 0;
}
}  // namespace

extern "C" int reactioncmp(CChemicalReaction rxn, CChemicalReaction rxn2) {
  auto *rxnm = (ChemicalReaction *)rxn;
  auto *rxn2m = (ChemicalReaction *)rxn2;

  if (!rxnm) {
    if (!rxn2m) {
      return 0;
    }
    return -1;
  }
  if (!rxn2m) {
    return 1;
  }

  int res = rxnm->getNumReactantTemplates() - rxn2m->getNumReactantTemplates();
  if (res) {
    return res;
  }
  res = rxnm->getNumProductTemplates() - rxn2m->getNumProductTemplates();
  if (res) {
    return res;
  }
  if (!getIgnoreReactionAgents()) {
    res = rxnm->getNumAgentTemplates() - rxn2m->getNumAgentTemplates();
    if (res) {
      return res;
    }
  }

  MoleculeDescriptors *rxn_react =
      calcMolecularDescriptorsReaction(rxnm, Reactant);
  MoleculeDescriptors *rxn2_react =
      calcMolecularDescriptorsReaction(rxn2m, Reactant);
  res = compareMolDescriptors(*rxn_react, *rxn2_react);
  delete (rxn_react);
  delete (rxn2_react);
  if (res) {
    return res;
  }
  MoleculeDescriptors *rxn_product =
      calcMolecularDescriptorsReaction(rxnm, Product);
  MoleculeDescriptors *rxn2_product =
      calcMolecularDescriptorsReaction(rxn2m, Product);
  res = compareMolDescriptors(*rxn_product, *rxn2_product);
  delete (rxn_product);
  delete (rxn2_product);
  if (res) {
    return res;
  }
  if (!getIgnoreReactionAgents()) {
    MoleculeDescriptors *rxn_agent =
        calcMolecularDescriptorsReaction(rxnm, Agent);
    MoleculeDescriptors *rxn2_agent =
        calcMolecularDescriptorsReaction(rxn2m, Agent);
    res = compareMolDescriptors(*rxn_agent, *rxn2_agent);
    delete (rxn_agent);
    delete (rxn2_agent);
    if (res) {
      return res;
    }
  }

  RDKit::MatchVectType matchVect;
  if (hasReactionSubstructMatch(*rxnm, *rxn2m, (!getIgnoreReactionAgents()))) {
    return 0;
  }
  return -1;
}

extern "C" CSfp makeReactionDifferenceSFP(CChemicalReaction data, int size,
                                          int fpType) {
  auto *rxn = (ChemicalReaction *)data;
  SparseFP *res = nullptr;

  try {
    if (fpType > 3 || fpType < 1) {
      elog(ERROR, "makeReactionDifferenceSFP: Unknown Fingerprint type");
    }
    RDKit::ReactionFingerprintParams params;
    params.fpType = static_cast<FingerprintType>(fpType);
    params.fpSize = size;
    params.includeAgents = (!getIgnoreReactionAgents());
    params.agentWeight = getReactionDifferenceFPWeightAgents();
    params.nonAgentWeight = getReactionDifferenceFPWeightNonagents();
    res = (SparseFP *)RDKit::DifferenceFingerprintChemReaction(*rxn, params);
  } catch (...) {
    elog(ERROR, "makeReactionDifferenceSFP: Unknown exception");
  }
  return (CSfp)res;
}

extern "C" CBfp makeReactionBFP(CChemicalReaction data, int size, int fpType) {
  auto *rxn = (ChemicalReaction *)data;
  ExplicitBitVect *res = nullptr;

  try {
    if (fpType > 5 || fpType < 1) {
      elog(ERROR, "makeReactionBFP: Unknown Fingerprint type");
    }
    RDKit::ReactionFingerprintParams params;
    params.fpType = static_cast<FingerprintType>(fpType);
    params.fpSize = size;
    params.includeAgents = (!getIgnoreReactionAgents());
    params.bitRatioAgents = getReactionStructuralFPAgentBitRatio();
    res = (ExplicitBitVect *)RDKit::StructuralFingerprintChemReaction(*rxn,
                                                                      params);
  } catch (...) {
    elog(ERROR, "makeReactionBFP: Unknown exception");
  }

  if (res) {
    std::string *sres = new std::string(BitVectToBinaryText(*res));
    delete res;
    return (CBfp)sres;
  } else {
    return nullptr;
  }
}
