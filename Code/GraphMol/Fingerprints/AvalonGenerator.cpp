//
//  Copyright (C) 2026 RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//  This is a pure RDKit port of the fingerprinting algorithm from the Avalon
//  toolkit (SetFingerprintBits() / SetFingerprintCountsWithFocus() in
//  ssmatch.c). The feature enumeration and hashing follow the original code
//  as closely as possible, but molecule perception (aromaticity, rings,
//  hydrogen counts) is done with RDKit's own machinery.
//  Comments were mainly preserved from the original Avalon code.
//  Github Copilot was used to assist in the porting process.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/Fingerprints/AvalonGenerator.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <RDGeneral/Invariant.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace RDKit {
namespace AvalonFP {
namespace {

// ---- minimal stand-ins for the data structures used by the Avalon code ----
constexpr int NONE = 0;
constexpr int ZERO_COUNT = 1;
constexpr int SINGLE = 1;
constexpr int DOUBLE = 2;
constexpr int TRIPLE = 3;
constexpr int AROMATIC = 4;
constexpr int ANY_BOND = 8;
constexpr int SUB_ONE = 1;
constexpr int SUB_MORE = 6;
constexpr int SUB_AS_IS = -2;

constexpr int USE_RING_PATTERN = 0x000001;
constexpr int USE_RING_PATH = 0x000002;
constexpr int USE_ATOM_SYMBOL_PATH = 0x000004;
constexpr int USE_ATOM_CLASS_PATH = 0x000008;
constexpr int USE_ATOM_COUNT = 0x000010;
constexpr int USE_AUGMENTED_ATOM = 0x000020;
constexpr int USE_HCOUNT_PATH = 0x000040;
constexpr int USE_HCOUNT_CLASS_PATH = 0x000080;
constexpr int USE_HCOUNT_PAIR = 0x000100;
constexpr int USE_BOND_PATH = 0x000200;
constexpr int USE_AUGMENTED_BOND = 0x000400;
constexpr int USE_RING_SIZE_COUNTS = 0x000800;
constexpr int USE_DEGREE_PATH = 0x001000;
constexpr int USE_CLASS_SPIDERS = 0x002000;
constexpr int USE_FEATURE_PAIRS = 0x004000;
constexpr int USE_SCAFFOLD_IDS = 0x100000;
constexpr int USE_SCAFFOLD_COLORS = 0x200000;
constexpr int USE_SCAFFOLD_LINKS = 0x400000;
constexpr int USE_NON_SSS_BITS = 0xF00000;

struct AvalonState {
  explicit AvalonState(const ROMol &mol) {
    atoms.reserve(mol.getNumAtoms());
    atomColors.resize(mol.getNumAtoms());
    atomRingFlags.resize(mol.getNumAtoms());
    atomSubDescriptors.resize(mol.getNumAtoms());
    for (const auto atom : mol.atoms()) {
      atoms.push_back(atom);
    }

    bonds.reserve(mol.getNumBonds());
    bondTypes.resize(mol.getNumBonds());
    bondColors.resize(mol.getNumBonds());
    bondRingFlags.resize(mol.getNumBonds());
    for (const auto bond : mol.bonds()) {
      bonds.push_back(bond);
      if (bond->getIsAromatic() || bond->getBondType() == Bond::AROMATIC) {
        bondTypes[bond->getIdx()] = 4;
      } else if (bond->getBondType() == Bond::SINGLE) {
        bondTypes[bond->getIdx()] = 1;
      } else if (bond->getBondType() == Bond::DOUBLE) {
        bondTypes[bond->getIdx()] = 2;
      } else if (bond->getBondType() == Bond::TRIPLE) {
        bondTypes[bond->getIdx()] = 3;
      } else {
        bondTypes[bond->getIdx()] = ANY_BOND;
      }
    }
  }

  std::vector<const Atom *> atoms;
  std::vector<const Bond *> bonds;
  std::vector<int> atomColors;
  std::vector<int> atomRingFlags;
  std::vector<int> atomSubDescriptors;
  std::vector<int> bondTypes;
  std::vector<int> bondColors;
  std::vector<int> bondRingFlags;
};

int &atomColor(AvalonState &state, const Atom *atom) {
  return state.atomColors[atom->getIdx()];
}

int &atomRingFlags(AvalonState &state, const Atom *atom) {
  return state.atomRingFlags[atom->getIdx()];
}

int &atomSubDescriptor(AvalonState &state, const Atom *atom) {
  return state.atomSubDescriptors[atom->getIdx()];
}

int &bondColor(AvalonState &state, const Bond *bond) {
  return state.bondColors[bond->getIdx()];
}

int &bondType(AvalonState &state, const Bond *bond) {
  return state.bondTypes[bond->getIdx()];
}

int &bondRingFlags(AvalonState &state, const Bond *bond) {
  return state.bondRingFlags[bond->getIdx()];
}

int bondEndpoint(const AvalonState &state, unsigned int bondIdx,
                 unsigned int end) {
  const auto bond = state.bonds[bondIdx];
  return static_cast<int>(
      (end == 0 ? bond->getBeginAtomIdx() : bond->getEndAtomIdx()) + 1);
}

struct neighbourhood_t {
  int n_ligands = 0;
  std::vector<int> atoms;  // 0-based atom indices
  std::vector<int> bonds;  // 0-based bond indices
};

// ---- hashing ----
uint64_t next_hash(uint64_t hash, uint64_t data) {
  hash += data;
  hash += (uint64_t)(hash << (uint64_t)10);
  hash ^= (uint64_t)(hash >> (uint64_t)6);
  return hash;
}

uint64_t hash_position(uint64_t hash, int nslots) {
  hash += (uint64_t)(hash << (uint64_t)3);
  hash ^= (uint64_t)(hash >> (uint64_t)11);
  hash += (uint64_t)(hash << (uint64_t)15);
  return (hash % (uint64_t)nslots);
}

// true if symbol is one of the comma separated tokens in list
bool AtomSymbolMatch(const std::string &symbol,
                     const std::vector<const char *> &list) {
  const auto symbol_cstr = symbol.c_str();
  return std::ranges::find_if(list, [symbol_cstr](const auto a) {
           return std::strcmp(a, symbol_cstr) == 0;
         }) != list.end();
}

// ---- the Avalon algorithm ----
constexpr int RING_PATTERN_SEED = 11;
constexpr int RING_PATH_SEED = 13;
constexpr int ATOM_SYMBOL_PATH_SEED = 17;
constexpr int ATOM_CLASS_PATH_SEED = 23;
constexpr int ATOM_COUNT_SEED = 31;
constexpr int AUGMENTED_ATOM_SEED = 37;
constexpr int HCOUNT_PATH_SEED = 41;
constexpr int HCOUNT_CLASS_PATH_SEED = 43;
constexpr int HCOUNT_PAIR_SEED = 47;
constexpr int BOND_PATH_SEED = 53;
constexpr int AUGMENTED_BOND_SEED = 61;
constexpr int RING_SIZE_SEED = 67;
constexpr int DEGREE_PATH_SEED = 71;
constexpr int CLASS_SPIDER_SEED = 79;
constexpr int RING_CLOSURE_SEED = 101;
constexpr int NON_SSS_SEED = 179;
#define MIN(a, b) ((a) < (b) ? (a) : (b))

/* new macro to convert the current seed value into the 'incremented' one */
#define NEXT_SEED(seed, increment) next_hash(seed, increment)
#define ADD_BIT(counts, ncounts, seed) (counts[hash_position(seed, ncounts)]++)
#define ADD_BIT_COUNT(counts, ncounts, seed, count) \
  (counts[hash_position(seed, ncounts)] += count)
#define SET_BIT(bytes, nbytes, seed) \
  (bytes[((seed) / 8) % nbytes] |= 0xFF & (1 << ((seed) % 8)))

/* Flags to be used to control recursive processing */
constexpr int PROCESS_RING_CLOSURES = 0x0001;
constexpr int PROCESS_CHAINS = 0x0002;
constexpr int FORCED_HETERO_END = 0x0004;
constexpr int IGNORE_PATH_SYMBOL = 0x0008;
constexpr int IGNORE_TERM_SYMBOL = 0x0010;
constexpr int FORCED_RING_PATH = 0x0020;
constexpr int STOP_AT_HEAVY_ATOM = 0x0040;
constexpr int DEBUG_PATH = 0x0100;

constexpr int ANY_COLOR = 113;

constexpr int CSP3 = 19;
constexpr int HETERO = 23;
constexpr int GENERIC = -1;

constexpr int SPECIAL_RING = (0xFC & ~(1 << 6));

static void SetPathLengthFlags(const ROMol &mol, AvalonState &state,
                               std::vector<int> &touched_indices,
                               int start_index, int path_length,
                               int current_index, int max_size,
                               std::vector<std::vector<int>> &length_matrix,
                               const std::vector<neighbourhood_t> &nbp,
                               int exclude_atom)
/*
 * Recursively traces the neighbouring of an atom (start_index+1)
 * collecting the path_lengths as bit flags in length_matrix[][].
 */
{
  for (int i = 0; i < nbp[current_index].n_ligands; i++) {
    if (path_length + 1 > max_size) { /* don't go too far */
      continue;
    }
    const auto ai = nbp[current_index].atoms[i];
    if (ai + 1 == exclude_atom) {
      continue;
    }
    if (touched_indices[ai]) {
      continue; /* don't walk backwards */
    }
    if (atomColor(state, state.atoms[ai]) == 0) {
      continue;
    }
    touched_indices[ai] = 1; /* updating */
    length_matrix[start_index][ai] |= 1 << (path_length + 1);
    SetPathLengthFlags(mol, state, touched_indices, start_index,
                       path_length + 1, ai, max_size, length_matrix, nbp,
                       exclude_atom);
    touched_indices[ai] = 0; /* down-dating */
  }
}

static void SpecialNeighboursRec(
    const ROMol &mol, AvalonState &state, std::vector<int> &touched_indices,
    int path_length, int current_index, int max_size,
    /* count of sp3 carbons with >= 3 C neighbours */
    int csp3[],
    /* count of hetero atoms */
    int hetero[], const std::vector<neighbourhood_t> &nbp, int exclude_atom)
/*
 * Recursively traces the neighbouring of an atom collecting
 * counts of special atoms at certain graph distances.
 *
 * exclude_atom terminates neighbourhood search paths.
 */
{
  for (int i = 0; i < nbp[current_index].n_ligands; i++) {
    if (path_length + 1 > max_size) { /* don't go too far */
      continue;
    }
    const auto ai = nbp[current_index].atoms[i];
    if (ai + 1 == exclude_atom) {
      continue;
    }
    if (touched_indices[ai]) {
      continue; /* don't walk backwards */
    }
    if (atomColor(state, state.atoms[ai]) <= 0) {
      continue;
    }
    if (atomColor(state, state.atoms[ai]) == CSP3) {
      csp3[path_length]++;
    } else if (atomColor(state, state.atoms[ai]) == HETERO) {
      hetero[path_length]++;
    }
    /* only carbon-connected SPIDERS (beta atoms are always included) */
    if (atomColor(state, state.atoms[ai]) == HETERO && path_length > 1) {
      continue;
    }
    touched_indices[ai] = 1; /* updating */
    SpecialNeighboursRec(mol, state, touched_indices, path_length + 1, ai,
                         max_size, csp3, hetero, nbp, exclude_atom);
    touched_indices[ai] = 0; /* down-dating */
  }
}
int SetPathBitsRec(const ROMol &mol, AvalonState &state,
                   const std::vector<neighbourhood_t> &nbp, int *fp_counts,
                   int ncounts, uint64_t seed,
                   std::vector<int> &touched_indices, int nbonds, int minbonds,
                   int maxbonds, int sprout_index, int first_index,
                   int last_index, int flags, int exclude_atom)
/*
 * Recursively enumerates the paths through *mp. The next sprouting
 * step is done on the atom (sprout_index+1). seed represents the
 * hash_code processing value up to and including this atom. touched_indices[i]
 * is > 0 if atom (i+1) is already included in the parent path.
 *
 * Paths through exclude_atom are terminated.
 */
{
  int result = 0;
  auto old_seed = seed;
  if (nbonds > maxbonds) {
    return result;
  }
  for (int i = 0; i < nbp[sprout_index].n_ligands; i++) {
    const auto ai = nbp[sprout_index].atoms[i];
    if (ai == last_index) {
      continue;
    }
    if (ai + 1 == exclude_atom) {
      continue;
    }
    const auto bi = nbp[sprout_index].bonds[i];
    const auto ap = mol.getAtomWithIdx(ai);
    if ((flags & FORCED_RING_PATH) &&
        bondRingFlags(state, mol.getBondWithIdx(bi)) == 0) {
      continue;
    }
    if (atomColor(state, ap) >= 18 && atomColor(state, ap) != ANY_COLOR &&
        0 != (flags & STOP_AT_HEAVY_ATOM)) {
      continue;
    }
    if (touched_indices[ai] > 0) /* ring closure */
    {
      if (0 == (flags & PROCESS_RING_CLOSURES)) {
        continue;
      }
      if (touched_indices[ai] > 1) {
        continue;  // not just a plain ring
      }
      const auto bcolor = bondColor(state, mol.getBondWithIdx(bi));
      if (bcolor == 0) {
        continue;
      }
      const auto acolor = atomColor(state, mol.getAtomWithIdx(ai));
      if (acolor == 0 && 0 == (flags & IGNORE_PATH_SYMBOL)) {
        continue;
      }
      old_seed = seed;
      seed = NEXT_SEED(seed, bcolor * 16);
      if (0 == (flags & IGNORE_PATH_SYMBOL)) {
        seed = NEXT_SEED(seed, acolor);
      }
      seed = NEXT_SEED(
          seed, touched_indices[ai]);  // codes how far this closure points back
      ADD_BIT(fp_counts, ncounts, seed);
      seed = NEXT_SEED(seed,
                       RING_CLOSURE_SEED *
                           (nbonds - touched_indices[ai]));  // codes ring size
      ADD_BIT(fp_counts, ncounts, seed);
      result++;
      seed = old_seed;
    } else /* normal path */
    {
      const auto bcolor = bondColor(state, mol.getBondWithIdx(bi));
      if (bcolor == 0) {
        continue;
      }
      const auto acolor = atomColor(state, mol.getAtomWithIdx(ai));
      if (acolor == 0 && 0 == (flags & IGNORE_PATH_SYMBOL)) {
        continue;
      }
      /* saving and updating */
      touched_indices[ai] = nbonds + 1;
      old_seed = seed;
      seed = NEXT_SEED(seed, bcolor * 16);
      if (0 == (flags & IGNORE_PATH_SYMBOL)) {
        seed = NEXT_SEED(seed, acolor);
      }
      if (nbonds >= minbonds &&
          (acolor > 0 || 0 != (flags & IGNORE_TERM_SYMBOL)) &&
          0 != (flags & PROCESS_CHAINS)) {
        if (0 == (flags & FORCED_HETERO_END) ||
            0 != (flags & IGNORE_TERM_SYMBOL) || acolor != 6) {
          if (acolor > 1) /* don't hit hydrogens! */
          {
            if (0 != (flags & IGNORE_TERM_SYMBOL)) {
              ADD_BIT(fp_counts, ncounts, seed);
            } else {
              ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, 17 * acolor));
            }
            result++;
          }
        }
      }
      if (nbonds + 1 <= maxbonds) /* continue recursion */
      {
        result +=
            SetPathBitsRec(mol, state, nbp, fp_counts, ncounts, seed,
                           touched_indices, nbonds + 1, minbonds, maxbonds, ai,
                           first_index, sprout_index, flags, exclude_atom);
      }

      /* restoring and down dating */
      touched_indices[ai] = 0;
      seed = old_seed;
    }
  }

  return result;
}

constexpr int HETERO_FLAG = 0x0100;
constexpr int RING_SUBST_FLAG = 0x0200;
constexpr int QUART_FLAG = 0x0400;
constexpr int CSP3_FLAG = 0x0800;
constexpr int RS_SPECIAL_FLAG = 0x1000;
constexpr int TYPE_MASK = 0x00FF;
constexpr int C_FLAG = 0x0001;
constexpr int O_FLAG = 0x0002;
constexpr int N_FLAG = 0x0003;
constexpr int S_FLAG = 0x0004;
constexpr int P_FLAG = 0x0005;
constexpr int X_FLAG = 0x0006;

int SetFeatureBits(const ROMol &mol, AvalonState &state, int *fp_counts,
                   int ncounts, int start_flags, int end_flags, int path_min,
                   int path_max, bool use_counts, bool use_atom_types,
                   const std::vector<std::vector<int>> &length_matrix,
                   uint64_t start_seed, int exclude_atom) {
  int result = 0;
  const auto nAtoms = static_cast<int>(mol.getNumAtoms());
  std::vector<int> counts(ncounts * 4, 0);
  for (int i = 0; i < nAtoms; i++) {
    if (i + 1 == exclude_atom) {
      continue;
    }
    const auto coli = atomColor(state, mol.getAtomWithIdx(i));
    if (0 == (coli & start_flags)) {
      continue;
    }
    if (use_atom_types && 0 == (coli & TYPE_MASK)) {
      continue;  // ignore generic atoms
    }
    const auto seed_i =
        use_atom_types ? NEXT_SEED(start_seed, coli & TYPE_MASK) : start_seed;
    for (int j = 0; j < nAtoms; j++) {
      if (j + 1 == exclude_atom) {
        continue;
      }
      const auto colj = atomColor(state, mol.getAtomWithIdx(j));
      if (0 == (colj & end_flags)) {
        continue;
      }
      if (use_atom_types && 0 == (colj & TYPE_MASK)) {
        continue;  // ignore generic atoms
      }
      const auto seed =
          use_atom_types ? NEXT_SEED(seed_i, colj & TYPE_MASK) : seed_i;
      for (int k = path_min; k <= path_max; k++) {
        if ((1 << k) & length_matrix[i][j]) {
          counts[(k * 19 + seed) % (ncounts * 4)]++;
          result++;
        }
      }
    }
  }
  /* Set the bits */
  for (int i = 0; i < ncounts * 4; i++) {
    if (counts[i] > 0) {
      ADD_BIT(fp_counts, ncounts, i);
      result++;
    }
    if (use_counts && counts[i] > 1) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i, 19), counts[i]);
      result++;
    }
    if (use_counts && counts[i] > 2) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i, 29), counts[i]);
      result++;
    }
    if (use_counts && counts[i] > 4) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i, 59), counts[i]);
      result++;
    }
    if (use_counts && counts[i] > 8) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i, 79), counts[i]);
      result++;
    }
  }
  return (result);
}
int countAtomCountFeatures(AvalonState &state, const int *H_count,
                           const int *atom_status, int *fp_counts, int nAtoms,
                           int ncounts, int exclude_atom, int nrare_atoms,
                           int ndouble, int naromatic, int nfusionb,
                           int which_bits) {
  int result = 0;
  auto ap = state.atoms.data();
  if (which_bits & USE_ATOM_COUNT) {
    constexpr int NCOUNT_HASH = 128;
    constexpr int NCOUNT_SEED_HASH = 128 * 128;
    std::array<int, NCOUNT_HASH> atom_type_count_hash{};
    std::array<int, NCOUNT_SEED_HASH> atom_type_count_seed_hash{};
    uint64_t seed = 0;
    int nbits = 0;
    int nringch2 = 0;
    int nfusionch = 0;
    int nspiro = 0;
    int hash = 0;
    /* Collect hashed counts of atom types with hydrogen counts */
    std::fill(atom_type_count_hash.begin(), atom_type_count_hash.end(), 0);
    std::fill(atom_type_count_seed_hash.begin(),
              atom_type_count_seed_hash.end(), 0);
    nringch2 = 0;
    nfusionch = 0;
    nspiro = 0;
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      if (i + 1 == exclude_atom) {
        continue;
      }
      if (atomColor(state, *ap) == 0) {
        continue;
      }
      if (atomColor(state, *ap) == 6) {
        /* count CH2 in ring */
        if (atom_status[i] > 0 && H_count[i + 1] >= 2) {
          nringch2++;
        }
        /* count ring fusion CH atoms */
        if (atom_status[i] > 2 && H_count[i + 1] >= 1) {
          nfusionch++;
        }
        if (atom_status[i] > 3) {
          nspiro++;
        }
        continue;  // carbon has a special count processing
      }
      hash = 0;
      seed = NEXT_SEED(317 * ATOM_COUNT_SEED, 507);
      seed = NEXT_SEED(seed, atomColor(state, *ap) + 17);
      atom_type_count_seed_hash[seed % NCOUNT_SEED_HASH]++;
      /* normal hetero only with hydrogen */
      if ((atomColor(state, *ap) == 7 || atomColor(state, *ap) == 8) &&
          H_count[i + 1] <= 0) {
        continue;
      }
      hash = hash * 7 + atomColor(state, *ap) + 13;
      atom_type_count_hash[hash % NCOUNT_HASH]++;
      seed = NEXT_SEED(seed, atomColor(state, *ap) + 13);
      atom_type_count_seed_hash[seed % NCOUNT_SEED_HASH]++;
      /* one more bit for rare types */
      if (atomColor(state, *ap) != 7 && atomColor(state, *ap) != 8) {
        hash = hash * 7 + atomColor(state, *ap) + 2 * 13;
        atom_type_count_hash[hash % NCOUNT_HASH]++;
        seed = NEXT_SEED(seed, atomColor(state, *ap) + 2 * 13);
        atom_type_count_seed_hash[seed % NCOUNT_SEED_HASH]++;
      }
    }
    nbits = 0;
    /* Now, we set the corresponding bits */

    if (ncounts <= 2048) {  // old atom count fingerprints
      for (int i = 0; i < NCOUNT_HASH; i++) {
        if (atom_type_count_hash[i] > 0) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 19, 3),
                        atom_type_count_hash[i]);
          nbits++;
        }
        if (atom_type_count_hash[i] > 1) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 23, 5),
                        atom_type_count_hash[i]);
          nbits++;
        }
        if (atom_type_count_hash[i] > 2) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 29, 7),
                        atom_type_count_hash[i]);
          nbits++;
        }
        if (atom_type_count_hash[i] > 4) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 31, 11),
                        atom_type_count_hash[i]);
          nbits++;
        }
        if (atom_type_count_hash[i] > 8) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 37, 13),
                        atom_type_count_hash[i]);
          nbits++;
        }
        if (atom_type_count_hash[i] > 16) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 41, 17),
                        atom_type_count_hash[i]);
          nbits++;
        }
      }
    } else {  // new atom count fingerprints
      for (int i = 0; i < NCOUNT_SEED_HASH; i++) {
        if (atom_type_count_seed_hash[i] > 0) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 19, 3),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
        if (atom_type_count_seed_hash[i] > 1) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 23, 5),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
        if (atom_type_count_seed_hash[i] > 2) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 29, 7),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
        if (atom_type_count_seed_hash[i] > 3) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 29, 71),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
        if (atom_type_count_seed_hash[i] > 4) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 31, 11),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
        if (atom_type_count_seed_hash[i] > 6) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 31, 113),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
        if (atom_type_count_seed_hash[i] > 8) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 37, 13),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
        if (atom_type_count_seed_hash[i] > 12) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 37, 133),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
        if (atom_type_count_seed_hash[i] > 16) {
          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i * 41, 17),
                        atom_type_count_seed_hash[i]);
          nbits++;
        }
      }
    }

    if (naromatic >= 10) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 341),
                    naromatic);
      nbits++;
    }
    if (naromatic >= 14) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 441),
                    naromatic);
      nbits++;
    }
    if (naromatic >= 18) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 541),
                    naromatic);
      nbits++;
    }
    if (naromatic >= 22) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 641),
                    naromatic);
      nbits++;
    }
    if (naromatic >= 26) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 741),
                    naromatic);
      nbits++;
    }

    if (ndouble >= 3) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 371),
                    naromatic);
      nbits++;
    }
    if (ndouble >= 5) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 471),
                    naromatic);
      nbits++;
    }
    if (ndouble >= 8) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 571),
                    naromatic);
      nbits++;
    }

    if (nringch2 >= 6) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 411),
                    nringch2);
      nbits++;
    }
    if (nringch2 >= 12) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 511),
                    nringch2);
      nbits++;
    }
    if (nringch2 >= 22) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 611),
                    nringch2);
      nbits++;
    }

    if (nfusionch >= 2) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 421),
                    nfusionch);
      nbits++;
    }
    if (nfusionch >= 4) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 521),
                    nfusionch);
      nbits++;
    }

    // Note spiro atoms
    if (nspiro >= 1) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 322),
                    nspiro);
      nbits++;
    }
    if (nspiro >= 1) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 422),
                    nspiro);
      nbits++;
    }
    if (nspiro >= 2) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 522),
                    nspiro);
      nbits++;
    }
    if (nspiro >= 2) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 622),
                    nspiro);
      nbits++;
    }

    if (nfusionb >= 1) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 425),
                    nfusionb);
      nbits++;
    }
    if (nfusionb >= 2) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 525),
                    nfusionb);
      nbits++;
    }
    if (nfusionb >= 3) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 625),
                    nfusionb);
      nbits++;
    }
    if (nfusionb >= 5) {
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 825),
                    nfusionb);
      nbits++;
    }
    if (nrare_atoms > 0) {
      seed = NEXT_SEED(317 * ATOM_COUNT_SEED, 5);
      ADD_BIT_COUNT(fp_counts, ncounts, seed, nrare_atoms);
      nbits++;
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(seed, 101), nrare_atoms);
      nbits++;
      ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(seed, 103), nrare_atoms);
      nbits++;
      if (nrare_atoms > 1) {
        ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(seed, 203), nrare_atoms);
        nbits++;
      }
      if (nrare_atoms > 3) {
        ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(seed, 211), nrare_atoms);
        nbits++;
      }
    }
    result += nbits;
  }

  return result;
}

int countAtomSymbolPathFeatures(const ROMol &mol, AvalonState &state,
                                const std::vector<neighbourhood_t> &nbp,
                                int *fp_counts, int ncounts,
                                std::vector<int> &touched_indices,
                                const std::vector<int> &degree,
                                const std::vector<int> &cdegree,
                                const int *atom_status, int which_bits,
                                int as_query, int exclude_atom) {
  const auto nAtoms = static_cast<int>(mol.getNumAtoms());
  int result = 0;
  uint64_t seed = 0;
  uint64_t old_seed = 0;
  auto ap = state.atoms.data();
  auto bp = state.bonds.data();
  if (which_bits & USE_ATOM_SYMBOL_PATH) {
    seed = ATOM_SYMBOL_PATH_SEED;
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      if (i + 1 == exclude_atom) {
        continue;
      }
      if (atomColor(state, *ap) <= 0) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, atomColor(state, *ap));
      /* Ignore common atom types */
      if (atomColor(state, *ap) != 6 && atomColor(state, *ap) != 7 &&
          atomColor(state, *ap) != 8) {
        ADD_BIT(fp_counts, ncounts, seed);
        result++;
      }
      /* fingerprint two atom pairs if substitution count is defined */
      if (1) {  // [TODO] 1
        if (atomColor(state, *ap) != 6 &&
            (!as_query || atomSubDescriptor(state, *ap) == SUB_AS_IS ||
             (atomSubDescriptor(state, *ap) != NONE &&
              atomSubDescriptor(state, *ap) != SUB_MORE &&
              atomSubDescriptor(state, *ap) == degree[i] + SUB_ONE - 1))) {
          result += SetPathBitsRec(
              mol, state, nbp, fp_counts, ncounts, seed + 12347,
              touched_indices, 1, 1, 1, /* path length 1 to 1 */
              i, 0, -1, FORCED_HETERO_END | PROCESS_CHAINS, exclude_atom);
        }
      }
      /* 2 bond paths for not very common starts */
      if (1) {  // [TODO] 2
        if (atomColor(state, *ap) >= 10 ||
            (atomColor(state, *ap) == 7 && degree[i] > 0) ||
            (atomColor(state, *ap) == 8 && degree[i] > 1) ||
            (atomColor(state, *ap) == 6 && degree[i] > 2 &&
             atom_status[i] > 0)) {
          result += SetPathBitsRec(
              mol, state, nbp, fp_counts, ncounts, seed, touched_indices, 1, 1,
              2, /* path length 1 to 2 */
              i, 0, -1,
              STOP_AT_HEAVY_ATOM |  // don't cross very heavy atoms
                  FORCED_HETERO_END | PROCESS_CHAINS,
              exclude_atom);
        }
      }
      /* Add more paths starting at special atoms */
      if (1) {  // [TODO] 3
        if (atomColor(state, *ap) > 6) {
          result += SetPathBitsRec(
              mol, state, nbp, fp_counts, ncounts,
              217 * atomColor(state, *ap) + seed, touched_indices, 1, 3,
              4, /* path length 3 to 4 */
              i, 0, -1,
              IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES |
                  STOP_AT_HEAVY_ATOM |  // don't cross very heavy atoms
                  PROCESS_CHAINS,
              exclude_atom);

          if (atomColor(state, *ap) > 10 &&
              atomColor(state, *ap) <=
                  18) {  // only third row of periodic table
            result += SetPathBitsRec(
                mol, state, nbp, fp_counts, ncounts, 17 + seed, touched_indices,
                1, 5, 7, i, 0, -1,
                FORCED_HETERO_END | IGNORE_PATH_SYMBOL |
                    STOP_AT_HEAVY_ATOM |  // don't cross very heavy atoms
                    // PROCESS_RING_CLOSURES |      // rings containing metal
                    // cause too many bits
                    PROCESS_CHAINS,
                exclude_atom);
          }
        }
      }
      if (1) {  // [TODO] 4
        if ((atomColor(state, *ap) == 7 ||
             atomColor(state, *ap) ==
                 8) &&  // only do this for common hetero elements
            cdegree[i] > 2) {
          seed = NEXT_SEED(seed, atomColor(state, *ap) * 23);
          result += SetPathBitsRec(
              mol, state, nbp, fp_counts, ncounts, seed, touched_indices, 1, 2,
              6, /* path length 1 to 4 */
              i, 0, -1,
              STOP_AT_HEAVY_ATOM |  // don't cross very heavy atoms
                  IGNORE_TERM_SYMBOL |
                  // IGNORE_PATH_SYMBOL |
                  PROCESS_CHAINS,
              exclude_atom);
        }
      }

      seed = old_seed;
      /* Add another set of bits for paths starting at more crowded atoms */
      if (1) {  // [TODO] 5
        if (degree[i] >= 4 &&
            atomColor(state, *ap) != 5 &&  // don't do it for boron
            atomColor(state, *ap) < 18)    // don't do it for transition metals
        {
          seed = NEXT_SEED(seed, atomColor(state, *ap));
          if (1) {
            result += SetPathBitsRec(
                mol, state, nbp, fp_counts, ncounts, NEXT_SEED(seed, 3 * 107),
                touched_indices, 1, 2, 4, /* path length 2 to 4 */
                i, 0, -1,
                STOP_AT_HEAVY_ATOM |     // don't cross very heavy atoms
                    IGNORE_TERM_SYMBOL | /* gedeck.mol */
                    PROCESS_CHAINS,
                exclude_atom);
          }
        }
      }

      seed = old_seed;
      /* Add another set of bits for paths starting at spiro atoms */
      if (1) {  // [TODO] 6
        if (atom_status[i] >= 4 &&
            atomColor(state, *ap) != 5 &&  // don't do it for boron
            atomColor(state, *ap) < 18)    // don't do it for transition metals
        {
          seed = NEXT_SEED(seed, atomColor(state, *ap) + 55);
          // if (1)
          result += SetPathBitsRec(
              mol, state, nbp, fp_counts, ncounts, NEXT_SEED(seed, 3 * 109),
              touched_indices, 1, 2, 4, /* path length 1 to 4 */
              i, 0, -1,
              STOP_AT_HEAVY_ATOM |  // don't cross very heavy atoms
                  IGNORE_TERM_SYMBOL | IGNORE_PATH_SYMBOL |
                  PROCESS_RING_CLOSURES | PROCESS_CHAINS,
              exclude_atom);
        }
      }

      seed = old_seed;
      /* Add bits for paths starting with rare bond orders */
      if (1) {
        for (int j = 0; j < nbp[i].n_ligands; j++) {
          bp = &state.bonds[nbp[i].bonds[j]];
          const auto ai = nbp[i].atoms[j];
          if (ai + 1 == exclude_atom) {
            continue;
          }
          if (bondColor(state, *bp) == 0) {
            continue;
          }
          if (bondType(state, *bp) != DOUBLE &&
              bondType(state, *bp) != TRIPLE) {
            continue;
          }
          seed = old_seed;
          seed = NEXT_SEED(seed, atomColor(state, *ap));
          seed = NEXT_SEED(seed, bondColor(state, *bp) * 613);
          touched_indices[ai] = 1; /* updating */
          result += SetPathBitsRec(
              mol, state, nbp, fp_counts, ncounts, seed, touched_indices, 2, 5,
              5, ai, 0, i,
              STOP_AT_HEAVY_ATOM |  // don't cross very heavy atoms
                  IGNORE_PATH_SYMBOL |
                  // IGNORE_TERM_SYMBOL |
                  PROCESS_RING_CLOSURES | (0 * PROCESS_CHAINS),
              exclude_atom);
          touched_indices[ai] = 0; /* down-dating */
        }
      }

      seed = old_seed;
      touched_indices[i] = 0; /* down-dating */
    }
  }
  return result;
}

int countAugmentedAtomFeatures(const ROMol &mol, AvalonState &state,
                               const std::vector<neighbourhood_t> &nbp,
                               int *fp_counts, int ncounts,
                               const std::vector<int> &degree,
                               const std::vector<int> &nspecial,
                               const int *H_count, int which_bits, int as_query,
                               int exclude_atom, uint64_t &seed,
                               uint64_t &old_seed) {
  const auto nAtoms = static_cast<int>(mol.getNumAtoms());
  int result = 0;
  auto ap = state.atoms.data();
  if (which_bits & USE_AUGMENTED_ATOM) {
    /* Set bits for all triples of atoms connected to a common atom */
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      if (atomColor(state, *ap) <= 0) {
        continue;
      }
      if (i + 1 == exclude_atom) {
        continue;
      }

      /* Set bit for atoms with more than one double or a triple bond */
      if (nspecial[i] >= 2) {
        ADD_BIT(fp_counts, ncounts,
                NEXT_SEED(atomColor(state, *ap) * AUGMENTED_ATOM_SEED, 101));
        // ADD_BIT(fp_counts, ncounts,
        // NEXT_SEED(atomColor(state, *ap)*AUGMENTED_ATOM_SEED,301));
        result += 1;
      }
      old_seed = seed;
      // Add some bits for hydrogen counted or hetero central atoms with
      // hetero neighbours
      if ((H_count[i + 1] > 0 || atomColor(state, *ap) != 6) &&
          degree[i] >= 2) {
        for (int i1 = 0; i1 < nbp[i].n_ligands; i1++) {
          if (state.atomColors[nbp[i].atoms[i1]] == 0) {
            continue;
          }
          if (nbp[i].atoms[i1] + 1 == exclude_atom) {
            continue;
          }
          if (state.bondColors[nbp[i].bonds[i1]] == 0) {
            continue;
          }
          for (int i2 = i1 + 1; i2 < nbp[i].n_ligands; i2++) {
            seed = NEXT_SEED(AUGMENTED_ATOM_SEED, 97);
            seed = NEXT_SEED(seed, atomColor(state, *ap));
            if (state.atomColors[nbp[i].atoms[i2]] == 0) {
              continue;
            }
            if (nbp[i].atoms[i2] + 1 == exclude_atom) {
              continue;
            }
            if (state.bondColors[nbp[i].bonds[i2]] == 0) {
              continue;
            }
            if (state.atomColors[nbp[i].atoms[i1]] == 6 &&
                state.atomColors[nbp[i].atoms[i2]] == 6) {
              continue;
            }
            int sum = 0;
            sum += state.atomColors[nbp[i].atoms[i1]] *
                   state.bondColors[nbp[i].bonds[i1]];
            sum += state.atomColors[nbp[i].atoms[i2]] *
                   state.bondColors[nbp[i].bonds[i2]];
            int prod = 1;
            prod *= state.atomColors[nbp[i].atoms[i1]];
            prod &= 0xFFF;
            prod *= state.atomColors[nbp[i].atoms[i2]];
            prod &= 0xFFF;
            seed = NEXT_SEED(seed, sum);
            seed = NEXT_SEED(seed, prod);
            ADD_BIT(fp_counts, ncounts, seed);
            // seed = NEXT_SEED(seed, (sum*prod)&0xFFF);
            // ADD_BIT(fp_counts, ncounts, seed);
            result += 1;
            int nmulti = 0;
            if (state.bondColors[nbp[i].bonds[i1]] >= 2) {
              nmulti++;
            }
            if (state.bondColors[nbp[i].bonds[i2]] >= 2) {
              nmulti++;
            }
            // add bit for hetero-substituted unsaturated hetero atom if
            // degree>1 is well-defined => catch e.g. nitroso vs. nitro
            if (atomColor(state, *ap) != 6 && nmulti > 0 &&
                (!as_query || atomSubDescriptor(state, *ap) == SUB_AS_IS ||
                 (atomSubDescriptor(state, *ap) != NONE &&
                  atomSubDescriptor(state, *ap) != SUB_MORE &&
                  atomSubDescriptor(state, *ap) == degree[i] + SUB_ONE - 1))) {
              ADD_BIT(fp_counts, ncounts,
                      NEXT_SEED(seed, nmulti + 173 * degree[i]));
              result += 1;
            }
          }
        }
      }

      if (degree[i] <= 2) {
        continue;
      }
      for (int i1 = 0; i1 < nbp[i].n_ligands; i1++) {
        if (nbp[i].atoms[i1] + 1 == exclude_atom) {
          continue;
        }
        for (int i2 = i1 + 1; i2 < nbp[i].n_ligands; i2++) {
          if (nbp[i].atoms[i2] + 1 == exclude_atom) {
            continue;
          }
          for (int i3 = i2 + 1; i3 < nbp[i].n_ligands; i3++) {
            if (nbp[i].atoms[i3] + 1 == exclude_atom) {
              continue;
            }
            seed = AUGMENTED_ATOM_SEED;
            seed = NEXT_SEED(seed, atomColor(state, *ap));
            if (state.atomColors[nbp[i].atoms[i1]] == 0) {
              continue;
            }
            if (state.atomColors[nbp[i].atoms[i2]] == 0) {
              continue;
            }
            if (state.atomColors[nbp[i].atoms[i3]] == 0) {
              continue;
            }
            if (state.bondColors[nbp[i].bonds[i1]] == 0) {
              continue;
            }
            if (state.bondColors[nbp[i].bonds[i2]] == 0) {
              continue;
            }
            if (state.bondColors[nbp[i].bonds[i3]] == 0) {
              continue;
            }
            /* count hetero neighbours */
            int nqtmp = 0;
            if (state.atomColors[nbp[i].atoms[i1]] != 6) {
              nqtmp++;
            }
            if (state.atomColors[nbp[i].atoms[i2]] != 6) {
              nqtmp++;
            }
            if (state.atomColors[nbp[i].atoms[i3]] != 6) {
              nqtmp++;
            }

            /* make sure to add some bits for really odd ones */
            int nmulti = 0;
            if (state.bondColors[nbp[i].bonds[i1]] == 2) {
              nmulti++;
            }
            if (state.bondColors[nbp[i].bonds[i1]] == 3) {
              nmulti++;
            }
            if (state.bondColors[nbp[i].bonds[i2]] == 2) {
              nmulti++;
            }
            if (state.bondColors[nbp[i].bonds[i2]] == 3) {
              nmulti++;
            }
            if (state.bondColors[nbp[i].bonds[i3]] == 2) {
              nmulti++;
            }
            if (state.bondColors[nbp[i].bonds[i3]] == 3) {
              nmulti++;
            }

            int sum = 0;
            sum += state.atomColors[nbp[i].atoms[i1]] *
                   state.bondColors[nbp[i].bonds[i1]];
            sum += state.atomColors[nbp[i].atoms[i2]] *
                   state.bondColors[nbp[i].bonds[i2]];
            sum += state.atomColors[nbp[i].atoms[i3]] *
                   state.bondColors[nbp[i].bonds[i3]];
            int prod = 1;
            prod *= state.atomColors[nbp[i].atoms[i1]];
            prod &= 0xFFF;
            prod *= state.atomColors[nbp[i].atoms[i2]];
            prod &= 0xFFF;
            prod *= state.atomColors[nbp[i].atoms[i3]];
            prod &= 0xFFF;
            seed = NEXT_SEED(seed, sum);
            seed = NEXT_SEED(seed, prod);
            if ((nqtmp > 2 || nmulti >= 2) && atomColor(state, *ap) == 6) {
              ADD_BIT(fp_counts, ncounts, seed);
              ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, 73));
              result += 2;
            }
            seed = NEXT_SEED(seed, (sum * prod) & 0xFFF);
            ADD_BIT(fp_counts, ncounts, seed);
            result++;
            if (nmulti >= 2 ||
                atomColor(state, *ap) > 6)  // make sure R-NO2 is covered
            {
              ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, 53));
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

int countAugmentedBondFeatures(AvalonState &state,
                               const std::vector<neighbourhood_t> &nbp,
                               int *fp_counts, int ncounts,
                               const std::vector<int> &degree, int nBonds,
                               int which_bits, int exclude_atom,
                               uint64_t &seed) {
  int result = 0;
  auto bp = state.bonds.data();
  if (which_bits & USE_AUGMENTED_BOND) {
    /* Set bits for all bonds with both end-degrees > 2 */
    bp = state.bonds.data();
    for (int i = 0; i < nBonds; i++, bp++) {
      if (bondEndpoint(state, (*bp)->getIdx(), 0) == exclude_atom) {
        continue;
      }
      if (bondEndpoint(state, (*bp)->getIdx(), 1) == exclude_atom) {
        continue;
      }
      const auto ai1 = bondEndpoint(state, (*bp)->getIdx(), 0) - 1;
      const auto ai2 = bondEndpoint(state, (*bp)->getIdx(), 1) - 1;
      if (degree[ai1] <= 2) {
        continue;
      }
      if (degree[ai2] <= 2) {
        continue;
      }
      for (int i1 = 0; i1 < nbp[ai1].n_ligands; i1++) {
        for (int i2 = i1 + 1; i2 < nbp[ai1].n_ligands; i2++) {
          /* don't reuse current bond */
          if (nbp[ai1].bonds[i1] == i) {
            continue;
          }
          if (nbp[ai1].bonds[i2] == i) {
            continue;
          }
          if (nbp[ai1].atoms[i1] + 1 == exclude_atom) {
            continue;
          }
          if (nbp[ai1].atoms[i2] + 1 == exclude_atom) {
            continue;
          }
          /* don't use masked atoms/bonds */
          if (state.atomColors[nbp[ai1].atoms[i1]] == 0) {
            continue;
          }
          if (state.atomColors[nbp[ai1].atoms[i2]] == 0) {
            continue;
          }
          if (state.bondColors[nbp[ai1].bonds[i1]] == 0) {
            continue;
          }
          if (state.bondColors[nbp[ai1].bonds[i2]] == 0) {
            continue;
          }
          int sumi = 0;
          sumi += state.atomColors[nbp[ai1].atoms[i1]] *
                  state.bondColors[nbp[ai1].bonds[i1]];
          sumi += state.atomColors[nbp[ai1].atoms[i2]] *
                  state.bondColors[nbp[ai1].bonds[i2]];
          int prodi = 1;
          prodi *= state.atomColors[nbp[ai1].atoms[i1]] +
                   state.bondColors[nbp[ai1].bonds[i1]];
          prodi *= state.atomColors[nbp[ai1].atoms[i2]] +
                   state.bondColors[nbp[ai1].bonds[i2]];
          sumi &= 0x0FFF;
          prodi &= 0x0FFF;
          for (int j1 = 0; j1 < nbp[ai2].n_ligands; j1++) {
            for (int j2 = j1 + 1; j2 < nbp[ai2].n_ligands; j2++) {
              /* don't reuse current bond */
              if (nbp[ai2].bonds[j1] == i) {
                continue;
              }
              if (nbp[ai2].bonds[j2] == i) {
                continue;
              }
              if (nbp[ai2].atoms[j1] + 1 == exclude_atom) {
                continue;
              }
              if (nbp[ai2].atoms[j2] + 1 == exclude_atom) {
                continue;
              }
              /* don't use masked atoms/bonds */
              if (state.atomColors[nbp[ai2].atoms[j1]] == 0) {
                continue;
              }
              if (state.atomColors[nbp[ai2].atoms[j2]] == 0) {
                continue;
              }
              if (state.bondColors[nbp[ai2].bonds[j1]] == 0) {
                continue;
              }
              if (state.bondColors[nbp[ai2].bonds[j2]] == 0) {
                continue;
              }
              int sumj = 0;
              sumj += state.atomColors[nbp[ai2].atoms[j1]] *
                      state.bondColors[nbp[ai2].bonds[j1]];
              sumj += state.atomColors[nbp[ai2].atoms[j2]] *
                      state.bondColors[nbp[ai2].bonds[j2]];
              int prodj = 1;
              prodj *= state.atomColors[nbp[ai2].atoms[j1]] +
                       state.bondColors[nbp[ai2].bonds[j1]];
              prodj *= state.atomColors[nbp[ai2].atoms[j2]] +
                       state.bondColors[nbp[ai2].bonds[j2]];
              sumj &= 0x0FFF;
              prodj &= 0x0FFF;
              seed = AUGMENTED_BOND_SEED;
              seed = NEXT_SEED(seed, bondColor(state, *bp));
              seed = NEXT_SEED(seed, prodi + prodj);
              seed = NEXT_SEED(seed, sumi * sumj);
              ADD_BIT(fp_counts, ncounts, seed);
              result++;
            }
          }
        }
      }
    }
  }
  return result;
}

int countHydrogenPairFeatures(AvalonState &state, int *fp_counts, int ncounts,
                              const int *H_count,
                              const std::vector<int> &unsaturated, int nBonds,
                              int which_bits, int exclude_atom,
                              uint64_t &seed) {
  int result = 0;
  auto bp = state.bonds.data();
  if (1 * which_bits & USE_HCOUNT_PAIR) {
    /* generate bits for hydrogen counted described bonds */
    bp = state.bonds.data();
    for (int i = 0; i < nBonds; i++, bp++) {
      if (H_count[bondEndpoint(state, (*bp)->getIdx(), 0)] == 0 &&
          H_count[bondEndpoint(state, (*bp)->getIdx(), 1)] == 0) {
        continue;
      }
      if (bondEndpoint(state, (*bp)->getIdx(), 0) == exclude_atom) {
        continue;
      }
      if (bondEndpoint(state, (*bp)->getIdx(), 1) == exclude_atom) {
        continue;
      }
      /* Don't consider CC single bonds */
      if (state.atomColors[bondEndpoint(state, (*bp)->getIdx(), 0) - 1] == 6 &&
          state.atomColors[bondEndpoint(state, (*bp)->getIdx(), 1) - 1] == 6 &&
          bondColor(state, *bp) == 1) {
        continue;
      }
      /* Don't consider explicit AH bonds */
      if (state.atomColors[bondEndpoint(state, (*bp)->getIdx(), 0) - 1] == 0) {
        continue;
      }
      if (state.atomColors[bondEndpoint(state, (*bp)->getIdx(), 1) - 1] == 0) {
        continue;
      }
      for (int j1 = 0; j1 <= H_count[bondEndpoint(state, (*bp)->getIdx(), 0)];
           j1++) {
        for (int j2 = 0; j2 <= H_count[bondEndpoint(state, (*bp)->getIdx(), 1)];
             j2++) {
          if (j1 + j2 == 0) {
            continue; /* at least one hydrogen */
          }
          // unsaturation triggers bit like a hydrogen
          if (!unsaturated[bondEndpoint(state, (*bp)->getIdx(), 0) - 1] &&
              !unsaturated[bondEndpoint(state, (*bp)->getIdx(), 1) - 1] &&
              j1 * j2 == 0 && bondColor(state, *bp) <= 1) {
            continue;
          }
          if (j1 + j2 > 3) {
            continue; /* at most 3 hydrogens */
          }
          // seed = HCOUNT_PAIR_SEED + 53*(j1+j2) + 7*(j1+1)*(j2+1);
          seed = NEXT_SEED(HCOUNT_PAIR_SEED, 7 * (j1 + 1) * (j2 + 1));
          seed = NEXT_SEED(seed, 53 * (j1 + j2));
          seed = NEXT_SEED(seed, bondColor(state, *bp));
          seed = NEXT_SEED(
              seed,
              state.atomColors[bondEndpoint(state, (*bp)->getIdx(), 0) - 1] +
                  state
                      .atomColors[bondEndpoint(state, (*bp)->getIdx(), 1) - 1]);
          seed = NEXT_SEED(
              seed,
              state.atomColors[bondEndpoint(state, (*bp)->getIdx(), 0) - 1] *
                  state
                      .atomColors[bondEndpoint(state, (*bp)->getIdx(), 1) - 1]);
          ADD_BIT(fp_counts, ncounts, seed);
          result++;
          if (bondColor(state, *bp) > 1) {
            seed = NEXT_SEED(seed, 83);
            ADD_BIT(fp_counts, ncounts, seed);
            result++;
            /* Add more bits for really special pairs */
            if (j1 + j2 == 1 &&
                (state.atomColors[bondEndpoint(state, (*bp)->getIdx(), 0) -
                                  1] != 6 ||
                 state.atomColors[bondEndpoint(state, (*bp)->getIdx(), 1) -
                                  1] != 6) &&
                bondColor(state, *bp) <= 3) {
              seed = NEXT_SEED(seed, 91);
              ADD_BIT(fp_counts, ncounts, seed);
              seed = NEXT_SEED(seed, 97);
              ADD_BIT(fp_counts, ncounts, seed);
              seed = NEXT_SEED(seed, 103);
              ADD_BIT(fp_counts, ncounts, seed);
              result += 3;
            }
          }
        }
      }
    }
  }
  return result;
}

int countHydrogenPathFeatures(const ROMol &mol, AvalonState &state,
                              const std::vector<neighbourhood_t> &nbp,
                              int *fp_counts, int ncounts,
                              std::vector<int> &touched_indices,
                              const std::vector<int> &degree,
                              const int *H_count, int which_bits,
                              int exclude_atom, uint64_t &seed,
                              uint64_t &old_seed) {
  const auto nAtoms = static_cast<int>(mol.getNumAtoms());
  int result = 0;
  auto ap = state.atoms.data();
  if (which_bits & USE_HCOUNT_PATH) {
    /* generate a short path for each atom that has a hydrogen */
    seed = HCOUNT_PATH_SEED;
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      if (atomColor(state, *ap) <= 0) {
        continue;
      }
      if (i + 1 == exclude_atom) {
        continue;
      }
      if (H_count[i + 1] == 0) {
        continue;
      }
      /* don't consider non methyl carbon atoms */
      if (atomColor(state, *ap) == 6 && H_count[i + 1] < 3 && degree[i] < 3) {
        continue;
      }

      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, atomColor(state, *ap));
      if (atomColor(state, *ap) != 6) {
        result += SetPathBitsRec(mol, state, nbp, fp_counts, ncounts, seed,
                                 touched_indices, 1, 2, 5,
                                 /* path length 2 to 4 */  // EVG
                                 i, 0, -1,
                                 IGNORE_PATH_SYMBOL |
                                     // IGNORE_TERM_SYMBOL |
                                     PROCESS_CHAINS,
                                 exclude_atom);
      } else {
        if (degree[i] > 2)  // tertiary hydrogen
        {
          if (0) {  // Class disabled to save bit density
            result += SetPathBitsRec(
                mol, state, nbp, fp_counts, ncounts, NEXT_SEED(seed, 101),
                touched_indices, 1, 2, 2, /* path length 2 to 3 */
                i, 0, -1,
                IGNORE_PATH_SYMBOL | IGNORE_TERM_SYMBOL | PROCESS_CHAINS,
                exclude_atom);
          }
          result += SetPathBitsRec(
              mol, state, nbp, fp_counts, ncounts, NEXT_SEED(seed, 103),
              touched_indices, 1, 2, 5, /* path length 2 to 3 */
              i, 0, -1,
              IGNORE_PATH_SYMBOL | FORCED_HETERO_END | /* i-Pr...Q */
                  PROCESS_CHAINS,
              exclude_atom);
        }
        if (H_count[i + 1] >= 3)  // methyl
        {
          if (0) {  // Class disabled to save bit density
            result += SetPathBitsRec(
                mol, state, nbp, fp_counts, ncounts, NEXT_SEED(seed, 1103),
                touched_indices, 1, 2, 3, /* path length 2 to 3 */
                i, 0, -1, IGNORE_PATH_SYMBOL | PROCESS_CHAINS, exclude_atom);
          }
          result += SetPathBitsRec(
              mol, state, nbp, fp_counts, ncounts, NEXT_SEED(seed, 1105),
              touched_indices, 1, 3, 6, /* path length 4 to 4 */
              i, 0, -1,
              IGNORE_PATH_SYMBOL | FORCED_HETERO_END | /* Me...Q */
                  PROCESS_CHAINS,
              exclude_atom);
        }
      }
      touched_indices[i] = 0; /* down-dating */
      seed = old_seed;
      if (atomColor(state, *ap) == 6) {
        continue;
      }

      if (H_count[i + 1] > 1) /* catch the difference between NH and NH2! */
      {
        touched_indices[i] = 1; /* updating */
        old_seed = seed;
        seed = NEXT_SEED(seed, 113);
        seed = NEXT_SEED(seed, atomColor(state, *ap));
        result +=
            SetPathBitsRec(mol, state, nbp, fp_counts, ncounts, seed,
                           touched_indices, 1, 1, 5, /* path length 1 to 5 */
                           i, 0, -1,
                           // FORCED_HETERO_END |
                           IGNORE_PATH_SYMBOL | PROCESS_CHAINS, exclude_atom);
        touched_indices[i] = 0; /* down-dating */
        seed = old_seed;
      }

      for (int j = 1; j < H_count[i + 1]; j++) {
        seed = NEXT_SEED(seed, 61 * j);
        seed = NEXT_SEED(seed, atomColor(state, *ap));
        ADD_BIT(fp_counts, ncounts, seed);
        result++;
      }
      seed = old_seed;
    }
  }
  return result;
}

int countRingSizeFeatures(AvalonState &state, int *fp_counts, int ncounts,
                          int nBonds, int which_bits, int exclude_atom) {
  int result = 0;
  uint64_t seed = 0;
  auto bp = state.bonds.data();
  if (which_bits & USE_RING_SIZE_COUNTS) {
    std::array<std::array<int, 15>, 15> rscounts{};
    for (int j = 3; j < 10; j++) /* loop through ring_sizes */
    {
      int nrbonds = 0;
      bp = state.bonds.data();
      for (int i = 0; i < nBonds; i++, bp++) {
        if (bondEndpoint(state, (*bp)->getIdx(), 0) != exclude_atom &&
            bondEndpoint(state, (*bp)->getIdx(), 1) != exclude_atom &&
            (bondRingFlags(state, *bp) & (1 << j))) {
          nrbonds++;
        }
      }
      seed = RING_SIZE_SEED;
      seed = NEXT_SEED(seed, j * 13);
      for (int i = 1; i < 100; i *= 2) {
        if (nrbonds >= j * i) {
          seed = NEXT_SEED(seed, i);
          /* don't set a bit if just one 5- or 6-ring */
          if ((j != 6 && j != 5) || nrbonds > j * i) {
            ADD_BIT_COUNT(fp_counts, ncounts, seed, nrbonds);
            result++;
          }
          if (j != 6 && j != 5) /* set one more bit for odd rings */
          {
            seed = NEXT_SEED(seed, 17);
            ADD_BIT_COUNT(fp_counts, ncounts, seed, nrbonds);
            result++;
          }
        } else {
          break;
        }
      }
    }
    /* Set bits for different ring sizes connected by a bond */
    for (int j = 0; j < 15; j++) {   /* loop through ring_sizes */
      for (int k = 0; k < 15; k++) { /* loop through ring_sizes */
        rscounts[j][k] = 0;
      }
    }
    /* count connections */
    bp = state.bonds.data();
    for (int i = 0; i < nBonds; i++, bp++) {
      if (bondEndpoint(state, (*bp)->getIdx(), 0) == exclude_atom) {
        continue;
      }
      if (bondEndpoint(state, (*bp)->getIdx(), 1) == exclude_atom) {
        continue;
      }
      for (int j = 3; j < 15; j++) { /* loop through ring_sizes */
        for (int k = 3; k < 15; k++) /* loop through ring_sizes */
        {
          if (j == k) {
            continue;
          }
          if ((state
                   .atomRingFlags[bondEndpoint(state, (*bp)->getIdx(), 0) - 1] &
               (1 << j)) &&
              (state
                   .atomRingFlags[bondEndpoint(state, (*bp)->getIdx(), 1) - 1] &
               (1 << k))) {
            rscounts[j][k]++;
            rscounts[k][j]++;
          }
        }
      }
    }
    /* set bits */
    for (int j = 3; j < 9; j++) { /* loop through not too large ring_sizes */
      for (int k = j + 1; k < 9;
           k++) /* loop through not too large ring_sizes */
      {
        if (rscounts[j][k] == 0) {
          continue;
        }
        seed = 2 * RING_SIZE_SEED;
        seed = NEXT_SEED(seed, (j + k) * 11);
        seed = NEXT_SEED(seed, j * k * 13);
        seed = NEXT_SEED(seed, 19);
        ADD_BIT(fp_counts, ncounts, seed);
        result++;
        /* more special links => more bits */
        if (rscounts[j][k] == 1 || j == 6 || k == 6) {
          continue;
        }
        seed = 2 * RING_SIZE_SEED + 23;
        seed = NEXT_SEED(seed, (j + k) * 11);
        seed = NEXT_SEED(seed, j * k * 13);
        seed = NEXT_SEED(seed, 19);
        ADD_BIT(fp_counts, ncounts, seed);
        result++;
      }
    }
  }
  return result;
}

int countScaffoldFeatures(AvalonState &state,
                          const std::vector<neighbourhood_t> &nbp,
                          int *fp_counts, int ncounts, const int *atom_status,
                          const int *bond_status,
                          const std::vector<int> &degree,
                          const std::vector<std::vector<int>> &length_matrix,
                          int nAtoms, int which_bits, int as_query,
                          int exclude_atom) {
  int result = 0;
  uint64_t seed = 0;
  auto ap = state.atoms.data();
  /* Collect bits that describe scaffolds. Those bits cannot (yet?) be used for
   * SSS screening */
  /* The method first collects a variant of the extended connectivity but
   * derived only        */
  /* on the ring bond graph. The bits are then set by combining atom features
   * with the        */
  /* corresponding scaffold extended connectivity of that atom. In addition, we
   * add           */
  /* bits to code connectivity of scaffolds by certain numbers of bonds using
   * the             */
  /* previously collected length_matrix. */
  if (!as_query && (which_bits & (USE_NON_SSS_BITS))) {
    std::vector<int> extcon(nAtoms, 0);
    std::vector<int> extcon2(nAtoms, 0);
    seed = NON_SSS_SEED;
    /* initialized extended connectivity */
    ap = state.atoms.data();
    for (int j = 0; j < nAtoms; j++, ap++) {
      if (atom_status[j] <= 0) {
        continue;
      }
      extcon[j] = atomRingFlags(state, *ap);
    }
    /* propagate extended connectivity to neighbours for a few cycles */
    for (int i = 0; i < 32; i++) {
      for (int j = 0; j < nAtoms; j++) {
        extcon2[j] = 0;
      }
      for (int j = 0; j < nAtoms; j++) {
        /* skip non-ring atoms */
        if (atom_status[j] <= 0) {
          continue;
        }
        extcon2[j] = atom_status[j] * 3 + (extcon[j] * 0xF);
        int sum = 0;
        int prod = 0;
        for (int jj = 0; jj < nbp[j].n_ligands; jj++) {
          // only propagate through ring bonds
          if (bond_status[nbp[j].bonds[jj]] <= 0) {
            continue;
          }
          sum += extcon[nbp[j].atoms[jj]];
          prod *= 0xFF & (1 + extcon[nbp[j].atoms[jj]]);
        }
        extcon2[j] = ((extcon2[j] + (sum * 191) + prod) << 8) +
                     ((extcon2[j] & 0xFF0000) >> 16);
        extcon2[j] &= 0xFFFFFF;
      }
      for (int j = 0; j < nAtoms; j++) {
        extcon[j] = extcon2[j];
      }
    }

    /* propagate smallest hash to all members of ring system */
    for (;;) {
      bool changed = false;
      for (int j = 0; j < nAtoms; j++) {
        /* skip non-ring atoms */
        if (atom_status[j] <= 0) {
          continue;
        }
        for (int jj = 0; jj < nbp[j].n_ligands; jj++) {
          // only propagate through ring bonds
          if (bond_status[nbp[j].bonds[jj]] <= 0) {
            continue;
          }
          if (extcon[j] < extcon[nbp[j].atoms[jj]]) {
            changed = true;
            extcon[nbp[j].atoms[jj]] = extcon[j];
          }
        }
      }
      if (!changed) {
        break;
      }
    }

    /* Now, use extcon to set bits */
    ap = state.atoms.data();
    for (int j = 0; j < nAtoms; j++, ap++) {
      if (j + 1 == exclude_atom) {
        continue;
      }
      /* skip non-ring atoms */
      if (atom_status[j] <= 0) {
        continue;
      }
      if (which_bits & USE_SCAFFOLD_IDS) {
        // ADD_BIT(fp_counts, ncounts, extcon[j]*1013 + seed);
        ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, extcon[j] * 10013));
        // ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, extcon[j]*40013));
        result += 1;
      }

      if (which_bits & USE_SCAFFOLD_LINKS) {
        /* needs to be a substitution point */
        if (degree[j] <= atom_status[j]) {
          continue;
        }
        for (int jj = 0; jj < nAtoms; jj++) {
          if (jj + 1 == exclude_atom) {
            continue;
          }
          /* skip non-ring atoms */
          if (atom_status[jj] <= 0) {
            continue;
          }
          if (j == jj) {
            continue;
          }
          /* needs to be a substitution point */
          if (degree[jj] <= atom_status[jj]) {
            continue;
          }
          int pathLength = 12;
          for (int k = 0; k < 12; k++) {
            if (length_matrix[j][jj] == (1 << k)) {
              pathLength = k;
              break;
            }
          }
          if (pathLength >= 12) {
            continue;  // multiple paths => no chain connection
          }
          if (pathLength >= 3) {
            continue;  // only short one counts
          }
          // bit for ring system
          ADD_BIT(fp_counts, ncounts,
                  NEXT_SEED(NEXT_SEED(seed, 3 * pathLength),
                            extcon[j] * 1013 + extcon[jj] * 2003));
          // bit for position
          ADD_BIT(fp_counts, ncounts,
                  NEXT_SEED(NEXT_SEED(seed, 5 * pathLength),
                            extcon2[j] * 2013 + extcon[jj] * 1003));
          ADD_BIT(fp_counts, ncounts,
                  NEXT_SEED(NEXT_SEED(seed, 5 * pathLength),
                            extcon[j] * 2013 + extcon2[jj] * 1003));
          ADD_BIT(fp_counts, ncounts,
                  NEXT_SEED(NEXT_SEED(seed, 7 * pathLength),
                            extcon2[j] * 3013 + extcon2[jj] * 3003));
          result += 4;
        }
      }
    }
    if (which_bits & USE_SCAFFOLD_COLORS) {
      ap = state.atoms.data();
      for (int j = 0; j < nAtoms; j++, ap++) {
        if (j + 1 == exclude_atom) {
          continue;
        }
        if (atom_status[j] <= 0) {
          continue;
        }
        int tmp1 = 0;
        if ((*ap)->getSymbol() == "C") {
          tmp1 = 101;
        } else if ((*ap)->getSymbol() == "O") {
          tmp1 = 301;
        } else if ((*ap)->getSymbol() == "N") {
          tmp1 = 401;
        } else if ((*ap)->getSymbol() == "S") {
          tmp1 = 601;
        } else if ((*ap)->getSymbol() == "P") {
          tmp1 = 701;
        } else if (AtomSymbolMatch((*ap)->getSymbol(),
                                   {"F", "Cl", "Br", "I", "At"})) {
          tmp1 = 901;
        }
        extcon[j] = atomRingFlags(state, *ap) + tmp1;
      }
      for (int i = 0; i < 32; i++) {
        for (int j = 0; j < nAtoms; j++) {
          extcon2[j] = 0;
        }
        for (int j = 0; j < nAtoms; j++) {
          /* skip non-ring atoms */
          if (atom_status[j] <= 0) {
            continue;
          }
          extcon2[j] = atom_status[j] * 3 + (extcon[j] * 0xF);
          int sum = 0;
          int prod = 0;
          for (int jj = 0; jj < nbp[j].n_ligands; jj++) {
            // only propagate through ring bonds
            if (bond_status[nbp[j].bonds[jj]] <= 0) {
              continue;
            }
            sum += extcon[nbp[j].atoms[jj]];
            prod *= 0xFF & (1 + extcon[nbp[j].atoms[jj]]);
          }
          extcon2[j] = ((extcon2[j] + (sum * 191) + prod) << 8) +
                       ((extcon2[j] & 0xFF0000) >> 16);
          extcon2[j] &= 0xFFFFFF;
        }
        for (int j = 0; j < nAtoms; j++) {
          extcon[j] = extcon2[j];
        }
      }
      /* propagate smallest hash to all members of ring system */
      for (;;) {
        bool changed = false;
        for (int j = 0; j < nAtoms; j++) {
          /* skip non-ring atoms */
          if (atom_status[j] <= 0) {
            continue;
          }
          for (int jj = 0; jj < nbp[j].n_ligands; jj++) {
            // only propagate through ring bonds
            if (bond_status[nbp[j].bonds[jj]] <= 0) {
              continue;
            }
            if (extcon[j] < extcon[nbp[j].atoms[jj]]) {
              changed = true;
              extcon[nbp[j].atoms[jj]] = extcon[j];
            }
          }
        }
        if (!changed) {
          break;
        }
      }
      for (int j = 0; j < nAtoms; j++) {
        if (j + 1 == exclude_atom) {
          continue;
        }
        ADD_BIT(fp_counts, ncounts,
                NEXT_SEED(NEXT_SEED(seed, 4), extcon[j] * 3013));
        // ADD_BIT(fp_counts, ncounts, NEXT_SEED(NEXT_SEED(seed,4),
        // extcon[j]*7013)); ADD_BIT(fp_counts, ncounts, extcon[j]*16013 +
        // 4*seed);
        result += 1;
      }
    }
  }
  return result;
}

int countFeaturePairFeatures(const ROMol &mol, AvalonState &state,
                             int *fp_counts, int ncounts, int nAtoms,
                             const int *atom_status,
                             const std::vector<int> &degree,
                             const std::vector<int> &cdegree,
                             const std::vector<std::vector<int>> &length_matrix,
                             int which_bits, int exclude_atom) {
  int result = 0;
  uint64_t seed = 0;
  auto ap = state.atoms.data();
  if (which_bits & USE_FEATURE_PAIRS) {
    /* set feature flags in atom colors */
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      int flags = 0;
      if ((*ap)->getSymbol() == "C") {
        flags = C_FLAG;
      } else if ((*ap)->getSymbol() == "O") {
        flags = O_FLAG;
      } else if ((*ap)->getSymbol() == "N") {
        flags = N_FLAG;
      } else if ((*ap)->getSymbol() == "S") {
        flags = S_FLAG;
      } else if ((*ap)->getSymbol() == "P") {
        flags = P_FLAG;
      } else if (AtomSymbolMatch((*ap)->getSymbol(),
                                 {"F", "Cl", "Br", "I", "At"})) {
        flags = X_FLAG;
      } else {
        flags = 0;
      }
      if (atomColor(state, *ap) == HETERO) {
        flags |= HETERO_FLAG;
      }
      if (cdegree[i] >= 3) {
        flags |= CSP3_FLAG;
      }
      if (degree[i] >= 4) {
        flags |= QUART_FLAG;
      }
      if (atom_status[i] > 0 && degree[i] >= 3) {
        flags |= RING_SUBST_FLAG;
        if (0 != (atomRingFlags(state, *ap) & SPECIAL_RING)) {
          flags |= RS_SPECIAL_FLAG;
        }
      }
      if (i + 1 == exclude_atom) {
        flags = 0;
      }
      atomColor(state, *ap) = flags;
    }
    /* collect bits for selected feature pairs */
    if (0) {  // Class disabled to save bit density
      result += SetFeatureBits(mol, state, fp_counts, ncounts, CSP3_FLAG,
                               HETERO_FLAG, /* to hetero or ring subst */
                               2, 3,        /* with path length 1 to 9 */
                               false,       /* don't use count */
                               true,        /* use atom type flags */
                               length_matrix, 1237, /* seed = 1237 */
                               exclude_atom);
    }
    if (0) {  // Class disabled to save bit density
      result += SetFeatureBits(mol, state, fp_counts, ncounts,
                               HETERO_FLAG, /* from ring substitution */
                               HETERO_FLAG, /* to hetero or ring subst */
                               1, 12,       /* with path length 1 to 10 */
                               true,        /* use count */
                               true,        /* use atom type flags */
                               length_matrix, 1237, /* seed = 1237 */
                               exclude_atom);
    }
    if (1) {
      result += SetFeatureBits(mol, state, fp_counts, ncounts,
                               RING_SUBST_FLAG, /* from ring substitution */
                               RING_SUBST_FLAG, /* to ring substitution */
                               5, 7,            /* with path length 1 to 12 */
                               false,           /* don't use count */
                               true,            /* use atom type flags */
                               length_matrix, 2237, /* seed = 2237 */
                               exclude_atom);
    }
    if (0) {  // Class disabled to save bit density
      result +=
          SetFeatureBits(mol, state, fp_counts, ncounts,
                         RS_SPECIAL_FLAG, /* from special ring substitution */
                         HETERO_FLAG,     /* to hetero atom */
                         2, 4,            /* with path length 1 to 5 */
                         true,            /* don't use count */
                         false,           /* use atom type flags */
                         length_matrix, 3237, /* seed = 3237 */
                         exclude_atom);
    }
    if (1) {
      result += SetFeatureBits(mol, state, fp_counts, ncounts,
                               QUART_FLAG,  /* from quartenary atom */
                               HETERO_FLAG, /* to hetero or ring subst */
                               1, 8,        /* with path length 1 to 8 */
                               false,       /* don't use count */
                               true,        /* use atom type flags */
                               length_matrix, 4237, /* seed = 4237 */
                               exclude_atom);
    }
    if (1) {
      result += SetFeatureBits(
          mol, state, fp_counts, ncounts, QUART_FLAG, /* from quartenary atom */
          RING_SUBST_FLAG, 1, 6, /* with path length 1 to 8 */
          false,                 /* don't use count */
          true,                  /* use atom type flags */
          length_matrix, 5237,   /* seed = 4237 */
          exclude_atom);
    }
    if (1) {
      result += SetFeatureBits(mol, state, fp_counts, ncounts,
                               X_FLAG,          /* from halogen atom */
                               CSP3_FLAG, 1, 1, /* with path length 1 to 1 */
                               false,           /* don't use count */
                               true,            /* use atom type flags */
                               length_matrix, 15237, /* seed = 15237 */
                               exclude_atom);
    }
    if (0) {  // Class disabled to save bit density
      result +=
          SetFeatureBits(mol, state, fp_counts, ncounts,
                         RS_SPECIAL_FLAG, /* from special ring substitution */
                         RING_SUBST_FLAG, /* to hetero atom */
                         2, 4,            /* with path length 1 to 5 */
                         true,            /* don't use count */
                         false,           /* use atom type flags */
                         length_matrix, 6237, /* seed = 6237 */
                         exclude_atom);
    }
    if (0) {  // Class disabled to save bit density
      result += SetFeatureBits(mol, state, fp_counts, ncounts,
                               HETERO_FLAG,     /* from hetero */
                               RING_SUBST_FLAG, /* to ring substitution */
                               1, 6,            /* with path length 1 to 8 */
                               false,           /* don't use count */
                               true,            /* use atom type flags */
                               length_matrix, 7237, /* seed = 4237 */
                               exclude_atom);
    }

    /* Set bits for ring-subst/ring-subst/hetero triples */
    if (0)  // too many spurious bits
    {
      auto ap1 = state.atoms.data();
      for (int i1 = 0; i1 < nAtoms; i1++, ap1++) {
        if (i1 + 1 == exclude_atom) {
          continue;
        }
        if (0 == (atomColor(state, *ap1) & RING_SUBST_FLAG)) {
          continue;
        }
        /* first atom must be in a non-sixmembered ring */
        if (!(atomRingFlags(state, *ap1) & SPECIAL_RING)) {
          continue;
        }
        auto ap2 = state.atoms.data();
        for (int i2 = 0; i2 < nAtoms; i2++, ap2++) {
          if (i1 == i2) {
            continue;
          }
          if (i2 + 1 == exclude_atom) {
            continue;
          }
          if (0 == (atomColor(state, *ap2) & RING_SUBST_FLAG)) {
            continue;
          }
          auto ap3 = state.atoms.data();
          for (int i3 = 0; i3 < nAtoms; i3++, ap3++) {
            if (i1 == i3) {
              continue;
            }
            if (i2 == i3) {
              continue;
            }
            if (i3 + 1 == exclude_atom) {
              continue;
            }
            if (0 == (atomColor(state, *ap3) & HETERO_FLAG)) {
              continue;
            }
            for (int j = 2; j <= 6; j++) {
              if (0 == (length_matrix[i1][i2] & (1 << j))) {
                continue;
              }
              for (int j1 = 2; j1 <= 5; j1++) {
                if (0 == (length_matrix[i2][i3] & (1 << j1))) {
                  continue;
                }
                for (int j2 = 2; j2 <= 5; j2++) {
                  if (0 == (length_matrix[i3][i1] & (1 << j2))) {
                    continue;
                  }
                  // Make sure we have a real triangle
                  if (j + j1 == j2) {
                    continue;
                  }
                  if (j + j2 == j1) {
                    continue;
                  }
                  if (j1 + j1 == j) {
                    continue;
                  }
                  if (j + j1 + j2 > 11) {
                    continue;  // spider too large
                  }
                  if (j + j1 + j2 < 8) {
                    continue;  // spider too small
                  }
                  seed = CLASS_SPIDER_SEED;
                  seed = NEXT_SEED(seed, j + j1 + j2);
                  seed = NEXT_SEED(seed, j * j1 * j2);
                  // distinguish ring sizes of second ring-subst
                  // for (k=3; k<15; k++)
                  for (int k = 3; k < 9; k++) {
                    if (!(atomRingFlags(state, *ap2) & (1 << k))) {
                      continue;
                    }
                    ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, k * 213));
                    result++;
                  }
                  // i1+1, atomRingFlags(state, *ap1), j,
                  // i2+1, atomRingFlags(state, *ap2), j1,
                  // i3+1, (*ap3)->getSymbol(), j2);
                }
              }
            }
          }
        }
      }
    }
  }
  return result;
}

int CountFingerprintPatterns(const ROMol &mol, AvalonState &state,
                             const std::vector<neighbourhood_t> &nbp,
                             int *H_count, int *atom_status, int *bond_status,
                             int *fp_counts, int ncounts, int which_bits,
                             int as_query, int exclude_atom) {
  const auto nAtoms = static_cast<int>(mol.getNumAtoms());
  const auto nBonds = static_cast<int>(mol.getNumBonds());
  std::vector<int> touched_indices(nAtoms, 0);
  std::vector<int> degree(nAtoms, 0);
  std::vector<int> cdegree(nAtoms, 0);
  std::vector<int> unsaturated(nAtoms, 0);
  std::vector<int> nspecial(nAtoms, 0);
  int nrare_atoms = 0;
  /* Set the color property to represent all different atom types */
  auto ap = state.atoms.data();
  for (int i = 0; i < nAtoms; i++, ap++) {
    unsaturated[i] = false;
    const auto symbol = (*ap)->getSymbol();
    atomColor(state, *ap) =
        symbol == "*" ? 0
                      : PeriodicTable::getTable()->getAtomicNumber(symbol);
    if (atomColor(state, *ap) <= 1) {
      atomColor(state, *ap) = 0; /* ignore hydrogens */
    }
    /* mark special atom types */
    if (atomColor(state, *ap) > 115) {
      atomColor(state, *ap) = -1;
    }
    if ((*ap)->getSymbol() == "A") {
      atomColor(state, *ap) = -1;
    }
    const auto is_rare =
        atomColor(state, *ap) > 0 &&
        !AtomSymbolMatch((*ap)->getSymbol(),
                         {"C", "H", "O", "N", "S", "P", "Cl", "F"});
    if (is_rare) {
      if (exclude_atom != i + 1 || exclude_atom <= 0) {
        nrare_atoms++;
      }
    }
  }

  int ndouble = 0;
  int naromatic = 0;
  int nfusionb = 0;
  /* Set the color property to represent the different bond type classes */
  auto bp = state.bonds.data();
  for (int i = 0; i < nBonds; i++, bp++) {
    if (bondType(state, *bp) == SINGLE) {
      bondColor(state, *bp) = 1;
    } else if (bondType(state, *bp) == DOUBLE) {
      bondColor(state, *bp) = 2;
      if (bondEndpoint(state, (*bp)->getIdx(), 0) != exclude_atom &&
          bondEndpoint(state, (*bp)->getIdx(), 1) != exclude_atom) {
        ndouble++;
      }
    } else if (bondType(state, *bp) == TRIPLE) {
      bondColor(state, *bp) = 3;
    } else if (bondType(state, *bp) == AROMATIC) {
      bondColor(state, *bp) = 4;
      if (bondEndpoint(state, (*bp)->getIdx(), 0) != exclude_atom &&
          bondEndpoint(state, (*bp)->getIdx(), 1) != exclude_atom) {
        naromatic++;
      }
    } else {
      bondColor(state, *bp) = 0;
    }
    if (bondColor(state, *bp) > 1) {
      unsaturated[bondEndpoint(state, (*bp)->getIdx(), 0) - 1] = true;
      unsaturated[bondEndpoint(state, (*bp)->getIdx(), 1) - 1] = true;
    }

    /* Count non-hydrogen degree */
    if (state.atomColors[bondEndpoint(state, (*bp)->getIdx(), 0) - 1] != 0 &&
        state.atomColors[bondEndpoint(state, (*bp)->getIdx(), 1) - 1] != 0) {
      degree[bondEndpoint(state, (*bp)->getIdx(), 0) - 1]++;
      degree[bondEndpoint(state, (*bp)->getIdx(), 1) - 1]++;
    }
    /* Count carbon degree */
    if (bondType(state, *bp) == DOUBLE) {
      nspecial[bondEndpoint(state, (*bp)->getIdx(), 0) - 1]++;
      nspecial[bondEndpoint(state, (*bp)->getIdx(), 1) - 1]++;
    } else if (bondType(state, *bp) == TRIPLE) {
      nspecial[bondEndpoint(state, (*bp)->getIdx(), 0) - 1] += 2;
      nspecial[bondEndpoint(state, (*bp)->getIdx(), 1) - 1] += 2;
    }
    if (bondType(state, *bp) != SINGLE) {
      continue;
    }
    if (state.atomColors[bondEndpoint(state, (*bp)->getIdx(), 0) - 1] != 0 &&
        state.atomColors[bondEndpoint(state, (*bp)->getIdx(), 1) - 1] == 6) {
      cdegree[bondEndpoint(state, (*bp)->getIdx(), 0) - 1]++;
    }
    if (state.atomColors[bondEndpoint(state, (*bp)->getIdx(), 1) - 1] != 0 &&
        state.atomColors[bondEndpoint(state, (*bp)->getIdx(), 0) - 1] == 6) {
      cdegree[bondEndpoint(state, (*bp)->getIdx(), 1) - 1]++;
    }

    if (exclude_atom <= 0 ||
        (bondEndpoint(state, (*bp)->getIdx(), 0) != exclude_atom &&
         bondEndpoint(state, (*bp)->getIdx(), 1) != exclude_atom)) {
      if (atom_status[bondEndpoint(state, (*bp)->getIdx(), 0) - 1] > 2 &&
          atom_status[bondEndpoint(state, (*bp)->getIdx(), 1) - 1] >
              2) {  // ring fusion
        nfusionb++;
      }
    }
  }
  /* ignore special atom types for further processing */
  ap = state.atoms.data();
  for (int i = 0; i < nAtoms; i++, ap++) {
    if (atomColor(state, *ap) < 0) {
      atomColor(state, *ap) = 0;
    }
  }

  int result = 0;
  if (which_bits & USE_ATOM_COUNT) {
    result += countAtomCountFeatures(state, H_count, atom_status, fp_counts,
                                     nAtoms, ncounts, exclude_atom, nrare_atoms,
                                     ndouble, naromatic, nfusionb, which_bits);
  }
  uint64_t seed = 0;
  uint64_t old_seed = 0;
  result += countAtomSymbolPathFeatures(
      mol, state, nbp, fp_counts, ncounts, touched_indices, degree, cdegree,
      atom_status, which_bits, as_query, exclude_atom);
  result += countAugmentedAtomFeatures(mol, state, nbp, fp_counts, ncounts,
                                       degree, nspecial, H_count, which_bits,
                                       as_query, exclude_atom, seed, old_seed);
  result += countAugmentedBondFeatures(state, nbp, fp_counts, ncounts, degree,
                                       nBonds, which_bits, exclude_atom, seed);
  result +=
      countHydrogenPairFeatures(state, fp_counts, ncounts, H_count, unsaturated,
                                nBonds, which_bits, exclude_atom, seed);
  result += countHydrogenPathFeatures(mol, state, nbp, fp_counts, ncounts,
                                      touched_indices, degree, H_count,
                                      which_bits, exclude_atom, seed, old_seed);
  /* Compute ring paths */
  ap = state.atoms.data();
  for (int i = 0; i < nAtoms; i++, ap++) {
    if (atom_status[i] <= 0) {
      atomColor(state, *ap) = 0;
    }
    if (i + 1 == exclude_atom) {
      atomColor(state, *ap) = 0;
    }
    if (atomColor(state, *ap) == 0) {
      continue;
    }
  }

  /* remove all bonds with only non-ring atoms from consideration */
  bp = state.bonds.data();
  for (int i = 0; i < nBonds; i++, bp++) {
    if (bond_status[i] <= 0 &&
        atom_status[bondEndpoint(state, (*bp)->getIdx(), 0) - 1] == 0 &&
        atom_status[bondEndpoint(state, (*bp)->getIdx(), 1) - 1] == 0) {
      bondColor(state, *bp) = 0;
    }
    if (bondEndpoint(state, (*bp)->getIdx(), 0) == exclude_atom) {
      bondColor(state, *bp) = 0;
    }
    if (bondEndpoint(state, (*bp)->getIdx(), 1) == exclude_atom) {
      bondColor(state, *bp) = 0;
    }
  }

  if (which_bits & USE_RING_PATH) {
    seed = RING_PATH_SEED;
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      if (atomColor(state, *ap) <= 0) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, atomColor(state, *ap));
      if (0) {  // Class disabled to save bit density
        result += SetPathBitsRec(mol, state, nbp, fp_counts, ncounts, seed,
                                 touched_indices, 1, 2, 3, i, 0, -1,
                                 PROCESS_CHAINS, exclude_atom);
      }

      if (atomColor(state, *ap) > 5 && atomColor(state, *ap) < 10 &&
          atom_status[i] > 2)  // only start at common light atoms
      {
        seed = NEXT_SEED(seed, 61);
        result +=
            SetPathBitsRec(mol, state, nbp, fp_counts, ncounts, seed,
                           touched_indices, 1, 3, 8, /* 3 to 8 */
                           i, 0, -1,
                           STOP_AT_HEAVY_ATOM |  // don't cross very heavy atoms
                               PROCESS_RING_CLOSURES,
                           exclude_atom);
      }
      seed = old_seed;

      touched_indices[i] = 0; /* down-dating */
    }
  }

  /* Set the color property to represent all different atom types */
  ap = state.atoms.data();
  for (int i = 0; i < nAtoms; i++, ap++) {
    const auto symbol = (*ap)->getSymbol();
    atomColor(state, *ap) =
        symbol == "*" ? 0
                      : PeriodicTable::getTable()->getAtomicNumber(symbol);
    if (atomColor(state, *ap) == 1 || i + 1 == exclude_atom) {
      atomColor(state, *ap) = 0; /* ignore hydrogens */
    } else {
      atomColor(state, *ap) = ANY_COLOR; /* treat all other atoms alike */
    }
  }

  /* Set the color property to represent the different bond type classes */
  bp = state.bonds.data();
  for (int i = 0; i < nBonds; i++, bp++) {
    if (bondType(state, *bp) == SINGLE) {
      bondColor(state, *bp) = 1;
    } else if (bondType(state, *bp) == DOUBLE) {
      bondColor(state, *bp) = 2;
    } else if (bondType(state, *bp) == TRIPLE) {
      bondColor(state, *bp) = 3;
    } else if (bondType(state, *bp) == AROMATIC) {
      bondColor(state, *bp) = 4;
    } else {
      bondColor(state, *bp) = 0;
    }
    if (bondEndpoint(state, (*bp)->getIdx(), 0) == exclude_atom) {
      bondColor(state, *bp) = 0;
    }
    if (bondEndpoint(state, (*bp)->getIdx(), 1) == exclude_atom) {
      bondColor(state, *bp) = 0;
    }
  }

  // Here, we have all non-trivial atoms mapped to ANY_COLOR while the
  // bond type is retained.
  if (which_bits & USE_BOND_PATH) {
    seed = BOND_PATH_SEED;
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      if (atomColor(state, *ap) <= 0) {
        continue;
      }
      if (i + 1 == exclude_atom) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      // start at branch node on ring
      if (degree[i] > 2 && atom_status[i] > 1) {
        result += SetPathBitsRec(mol, state, nbp, fp_counts, ncounts, seed,
                                 touched_indices, 1, 4, 4, /* was 4 to 4 */
                                 i, 0, -1, PROCESS_CHAINS, exclude_atom);
        if (atomRingFlags(state, *ap) & SPECIAL_RING) {
          result += SetPathBitsRec(
              mol, state, nbp, fp_counts, ncounts, NEXT_SEED(seed, 217),
              touched_indices, 1, 5, 5, i, 0, -1,
              IGNORE_PATH_SYMBOL | PROCESS_CHAINS, exclude_atom);
        }
      }
      /* Add other bits to catch poorly specified ring closures */
      seed = old_seed;
      seed = NEXT_SEED(seed, 11);
      result += SetPathBitsRec(
          mol, state, nbp, fp_counts, ncounts, seed, touched_indices, 1, 4,
          6, /* was 5 to 6 */
          i, 0, -1,
          // DEBUG_PATH |
          FORCED_RING_PATH | IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES,
          exclude_atom);
      /* Add bits for paths starting with rare bond orders */
      for (int j = 0; j < nbp[i].n_ligands; j++) {
        bp = &state.bonds[nbp[i].bonds[j]];
        const auto ai = nbp[i].atoms[j];
        if (ai + 1 == exclude_atom) {
          continue;
        }
        if (bondColor(state, *bp) == 0) {
          continue;
        }
        if (bondType(state, *bp) != DOUBLE && bondType(state, *bp) != TRIPLE) {
          continue;
        }
        if (atom_status[i] <= 0 && bondType(state, *bp) != TRIPLE) {
          continue;
        }
        seed = old_seed;
        seed = NEXT_SEED(seed, bondColor(state, *bp) * 413);
        touched_indices[ai] = 1; /* updating */
        result += SetPathBitsRec(
            mol, state, nbp, fp_counts, ncounts, seed, touched_indices, 2, 4, 5,
            ai, 0, i,
            IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES | PROCESS_CHAINS,
            exclude_atom);
        touched_indices[ai] = 0; /* down-dating */
      }
      seed = old_seed;
      touched_indices[i] = 0; /* down-dating */
    }
  }

  /* Set the color property to represent the different atom type classes */
  ap = state.atoms.data();
  for (int i = 0; i < nAtoms; i++, ap++) {
    if (i + 1 == exclude_atom) {
      atomColor(state, *ap) = 0;
      continue;
    }
    {
      if ((*ap)->getSymbol() == "H") {
        atomColor(state, *ap) = 0; /* ignore hydrogens */
      } else if ((*ap)->getSymbol() == "D") {
        atomColor(state, *ap) = 0; /* ignore hydrogens */
      } else if ((*ap)->getSymbol() == "T") {
        atomColor(state, *ap) = 0; /* ignore hydrogens */
      } else if ((*ap)->getSymbol() == "C") {
        atomColor(state, *ap) =
            6; /* carbon second row elements are one class */
      } else if ((*ap)->getSymbol() == "N") {
        atomColor(state, *ap) =
            8; /* nitrogen, oxigen, and sulfur are one class */
      } else if ((*ap)->getSymbol() == "O") {
        atomColor(state, *ap) =
            8; /* nitrogen, oxigen, and sulfur are one class */
      } else if ((*ap)->getSymbol() == "S") {
        atomColor(state, *ap) =
            8; /* nitrogen, oxigen, and sulfur are one class */
      } else if ((*ap)->getSymbol() == "Q") {
        atomColor(state, *ap) =
            8; /* nitrogen, oxigen, and sulfur are one class */
      } else if ((*ap)->getSymbol() == "A") {
        atomColor(state, *ap) = 0;
      } else {
        const auto symbol = (*ap)->getSymbol();
        const auto tmp =
            symbol == "*" ? 0
                          : PeriodicTable::getTable()->getAtomicNumber(symbol);
        if (1 < tmp && tmp < 115) {
          atomColor(state, *ap) = 8; /* non-carbon is the only second class */
        } else { /* This could be R atoms or other odd things */
          atomColor(state, *ap) = 0;
        }
      }
    }
    if (atomColor(state, *ap) > 115) {
      atomColor(state, *ap) = 0; /* ignore special atom types */
    }
  }

  // Here, we unify atom types to carbon/hetero distinction
  if (which_bits & USE_HCOUNT_CLASS_PATH) {
    /* generate a short path for each atom that has a hydrogen */
    seed = HCOUNT_CLASS_PATH_SEED;
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      if (atomColor(state, *ap) <= 0) {
        continue;
      }
      if (i + 1 == exclude_atom) {
        continue;
      }
      if (H_count[i + 1] == 0) {
        continue;
      }
      if (atomColor(state, *ap) == 6 && H_count[i + 1] < 2) {
        continue;
      }
      old_seed = seed;
      seed = NEXT_SEED(seed, atomColor(state, *ap));
      touched_indices[i] = 1; /* updating */
      if (atomColor(state, *ap) == 6) {
        result += SetPathBitsRec(
            mol, state, nbp, fp_counts, ncounts, seed, touched_indices, 1, 2,
            4, /* path length 1 to 3 */
            i, 0, -1, FORCED_HETERO_END | PROCESS_CHAINS, exclude_atom);
      } else {
        result += SetPathBitsRec(
            mol, state, nbp, fp_counts, ncounts, seed, touched_indices, 1, 2, 5,
            i, 0, -1, IGNORE_PATH_SYMBOL | FORCED_HETERO_END | PROCESS_CHAINS,
            exclude_atom);
      }
      seed = old_seed;
      touched_indices[i] = 0; /* down-dating */
    }
  }

  /* Do a small first pass with defined bond types */
  if (which_bits & USE_ATOM_CLASS_PATH) {
    // seed = ATOM_CLASS_PATH_SEED+117;
    seed = NEXT_SEED(ATOM_CLASS_PATH_SEED, 117);
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      if (atomColor(state, *ap) <= 0) {
        continue;
      }
      if (i + 1 == exclude_atom) {
        continue;
      }
      if (atomColor(state, *ap) == 6 && degree[i] < 3) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, atomColor(state, *ap));
      if (atomColor(state, *ap) == 6) {
        if (0) {  // Class disabled to save bit density
          result += SetPathBitsRec(
              mol, state, nbp, fp_counts, ncounts, seed, touched_indices, 1, 3,
              4, /* path length 3 to 4 */
              i, 0, -1, FORCED_HETERO_END | IGNORE_PATH_SYMBOL | PROCESS_CHAINS,
              exclude_atom);
        }
      } else {
        if (0) {  // Class disabled to save bit density
          result +=
              SetPathBitsRec(mol, state, nbp, fp_counts, ncounts, seed,
                             touched_indices, 1, 2, 3, /* path length 3 to 3 */
                             i, 0, -1,
                             // IGNORE_PATH_SYMBOL |
                             PROCESS_CHAINS, exclude_atom);
        }
        if (0) {  // Class disabled to save bit density
          result += SetPathBitsRec(
              mol, state, nbp, fp_counts, ncounts, seed, touched_indices, 1, 4,
              5, /* path length 4 to 5 */
              i, 0, -1, IGNORE_PATH_SYMBOL | FORCED_HETERO_END | PROCESS_CHAINS,
              exclude_atom);
        }
        result += SetPathBitsRec(
            mol, state, nbp, fp_counts, ncounts, seed, touched_indices, 1, 3,
            4, /* path length 4 to 7 */
            i, 0, -1, FORCED_RING_PATH | PROCESS_RING_CLOSURES | PROCESS_CHAINS,
            exclude_atom);
      }
      touched_indices[i] = 0;
      seed = old_seed;
    }
  }

  /* Set the color property to only a single class */
  bp = state.bonds.data();
  for (int i = 0; i < nBonds; i++, bp++) {
    if (SINGLE <= bondType(state, *bp) && bondType(state, *bp) <= ANY_BOND) {
      bondColor(state, *bp) = 5;
    } else {
      bondColor(state, *bp) = 0;
    }
    if (bondEndpoint(state, (*bp)->getIdx(), 0) == exclude_atom) {
      bondColor(state, *bp) = 0;
    }
    if (bondEndpoint(state, (*bp)->getIdx(), 1) == exclude_atom) {
      bondColor(state, *bp) = 0;
    }
  }

  // Here, we've unified atom types to carbon/hetero and made all bond types
  // identical
  if (which_bits & USE_ATOM_CLASS_PATH) {
    seed = ATOM_CLASS_PATH_SEED;
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      if (atomColor(state, *ap) <= 0) {
        continue;
      }
      if (i + 1 == exclude_atom) {
        continue;
      }
      if (atomColor(state, *ap) == 6) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, atomColor(state, *ap));
      // if (atomColor(state, *ap) == 6)
      {
        if (0 * degree[i] > 2) {  // Class disabled to save bit density
          result += SetPathBitsRec(
              mol, state, nbp, fp_counts, ncounts, seed, touched_indices, 1, 3,
              3, /* path length 3 to 3 */
              i, 0, -1, FORCED_HETERO_END | PROCESS_CHAINS, exclude_atom);
        }
      }
      // else
      {
        if (0) {  // Class disabled to save bit density
          result += SetPathBitsRec(
              mol, state, nbp, fp_counts, ncounts, seed, touched_indices, 1, 3,
              4, /* path length 3 to 4 */
              i, 0, -1, IGNORE_PATH_SYMBOL | PROCESS_CHAINS, exclude_atom);
        }
        seed = NEXT_SEED(seed, 23 + atomColor(state, *ap) * 19);
        result += SetPathBitsRec(
            mol, state, nbp, fp_counts, ncounts, seed, touched_indices, 1, 3,
            9, /* path length 3 to 9 */
            i, 0, -1, IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES, exclude_atom);
      }
      seed = old_seed;
      touched_indices[i] = 0; /* down-dating */
    }
    /* Q-Q and Q-C ring bond count */
    int qq_count = 0;
    int qc_count = 0;
    bp = state.bonds.data();
    for (int i = 0; i < nBonds; i++, bp++) {
      if (bond_status[i] == 0) {
        continue;
      }
      if (bondColor(state, *bp) == 0) {
        continue;
      }
      if (bondEndpoint(state, (*bp)->getIdx(), 0) == exclude_atom) {
        continue;
      }
      if (bondEndpoint(state, (*bp)->getIdx(), 1) == exclude_atom) {
        continue;
      }
      const auto ai1 =
          state.atomColors[bondEndpoint(state, (*bp)->getIdx(), 0) - 1];
      const auto ai2 =
          state.atomColors[bondEndpoint(state, (*bp)->getIdx(), 1) - 1];
      if (ai1 == 0 || ai2 == 0) {
        continue;
      }
      if (ai1 == 6 && ai2 == 6) {
        continue;
      }
      if (ai1 != 6 && ai2 != 6) {
        qq_count++;
        for (int j = 3; j < 9; j++) { /* set bits for not too large ring size */
          if (bondRingFlags(state, *bp) & (1 << j)) {
            ADD_BIT(fp_counts, ncounts,
                    NEXT_SEED(ATOM_CLASS_PATH_SEED * 17, j * 8));
            result++;
          }
        }
      } else {
        qc_count++;
        for (int j = 3; j < 9; j++) { /* set bits for not too large ring size */
          if (bondRingFlags(state, *bp) & (1 << j)) {
            ADD_BIT(fp_counts, ncounts,
                    NEXT_SEED(ATOM_CLASS_PATH_SEED * 19, j * 8));
            result++;
          }
        }
      }
    }
    seed = 2 * ATOM_CLASS_PATH_SEED + 3;
    for (int i = 1; i <= qq_count; i = (int)(1 + i * 1.5)) {
      seed = NEXT_SEED(seed, i * 153);
      ADD_BIT(fp_counts, ncounts, seed);
      result++;
      seed = NEXT_SEED(seed, 53);
      if (i <= 1) {
        ADD_BIT(fp_counts, ncounts, seed);
        result++;
      }
    }
    seed = 3 * ATOM_CLASS_PATH_SEED + 5;
    for (int i = 1; i <= MIN(qc_count, 2); i++) {
      seed = NEXT_SEED(seed, i * 157);
      ADD_BIT(fp_counts, ncounts, seed);
      result++;
    }
    for (int i = 3; i <= qc_count; i = (int)(i * 1.8)) {
      seed = NEXT_SEED(seed, i * 157);
      ADD_BIT(fp_counts, ncounts, seed);
      result++;
    }
  }

  /* Compute ring patters with at least one cycle */
  /* remove bonds from consideration that don't have at least one ring atom */
  bp = state.bonds.data();
  for (int i = 0; i < nBonds; i++, bp++) {
    if (bondEndpoint(state, (*bp)->getIdx(), 0) == exclude_atom) {
      continue;
    }
    if (bondEndpoint(state, (*bp)->getIdx(), 1) == exclude_atom) {
      continue;
    }
    bondColor(state, *bp) = 5;
    /* ignore non-ring bonds */
    if (atom_status[bondEndpoint(state, (*bp)->getIdx(), 0) - 1] == 0 &&
        atom_status[bondEndpoint(state, (*bp)->getIdx(), 1) - 1] == 0) {
      bondColor(state, *bp) = 0;
    } else {
      // may be redundant since atom types have already been unified
      ap = &state.atoms[bondEndpoint(state, (*bp)->getIdx(), 0) - 1];
      if (atomColor(state, *ap) != 0) {
        if (atomColor(state, *ap) != 6) {
          atomColor(state, *ap) = 8; /* non-carbons in one class */
        }
        if (atomRingFlags(state, *ap) == 0) {
          atomColor(state, *ap) = 0;
        }
      }
      // may be redundant since atom types have already been unified
      ap = &state.atoms[bondEndpoint(state, (*bp)->getIdx(), 1) - 1];
      if (atomColor(state, *ap) != 0) {
        if (atomColor(state, *ap) != 6) {
          atomColor(state, *ap) = 8; /* non-carbons in one class */
        }
        if (atomRingFlags(state, *ap) == 0) {
          atomColor(state, *ap) = 0;
        }
      }
    }
  }

  // Here, we have bond type ignored and atom types mapped to C and Q
  if (which_bits & USE_RING_PATTERN) {
    /* first process ring bond paths with atom classes */
    seed = RING_PATTERN_SEED;
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      if (atomColor(state, *ap) <= 0) {
        continue;
      }
      if (i + 1 == exclude_atom) {
        continue;
      }
      /* Don't process fragments starting at carbon in 6-ring only */
      if (atomColor(state, *ap) == 6 &&
          0 == (atomRingFlags(state, *ap) & SPECIAL_RING)) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, atomColor(state, *ap));
      result += SetPathBitsRec(mol, state, nbp, fp_counts, ncounts, seed,
                               touched_indices, 1, 3,
                               3, /* ring bond path size 3 to 3 */
                               i, 0, -1, PROCESS_CHAINS, exclude_atom);
      seed = old_seed;
      touched_indices[i] = 0; /* down-dating */
    }

    /* Now, we only include complete rings but ignore atom-type */
    /* 'A' atoms are now included nodes */
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      /* add 'A' atom to standard class */
      if ((*ap)->getSymbol() == "A") {
        atomColor(state, *ap) = 9;
      }
      if (i + 1 == exclude_atom) {
        atomColor(state, *ap) = 0;
      }
      if (atomColor(state, *ap) == 0) {
        continue;
      }
      atomColor(state, *ap) = 9; /* all ring atoms in same class */
    }
    bp = state.bonds.data();
    for (int i = 0; i < nBonds; i++, bp++) {
      if (bondEndpoint(state, (*bp)->getIdx(), 0) == exclude_atom) {
        bondColor(state, *bp) = 0;
      }
      if (bondEndpoint(state, (*bp)->getIdx(), 1) == exclude_atom) {
        bondColor(state, *bp) = 0;
      }
      if (bond_status[i] <= 0) {
        bondColor(state, *bp) = 0;
      }
    }
    // seed = RING_PATTERN_SEED+23;
    seed = NEXT_SEED(RING_PATTERN_SEED, 23);
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      if (atomColor(state, *ap) == 0) {
        continue;
      }
      if (i + 1 == exclude_atom) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, atomColor(state, *ap));
      if (0) {  // Class disabled to save bit density
        result += SetPathBitsRec(
            mol, state, nbp, fp_counts, ncounts, seed, touched_indices, 1, 4,
            17, /* ring size 4 to 17 */
            i, 0, -1, IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES, exclude_atom);
      }

      seed = old_seed;
      if (0) {                   // Class disabled to save bit density
        if (atom_status[i] > 2)  // start at ring fusion
        {
          seed = NEXT_SEED(seed, 61);
          result += SetPathBitsRec(
              mol, state, nbp, fp_counts, ncounts, seed, touched_indices, 1, 6,
              17, /* ring path size 6 to 17 */
              i, 0, -1, IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES,
              exclude_atom);
        }
      }

      seed = old_seed;
      if (0) {  // Class disabled to save bit density
        if (degree[i] > 2 && atom_status[i] >= 2)  // start at ring substituents
        {
          seed = NEXT_SEED(seed, 67);
          result += SetPathBitsRec(
              mol, state, nbp, fp_counts, ncounts, seed, touched_indices, 1, 6,
              17, /* ring path size 6 to 17 */
              i, 0, -1, IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES,
              exclude_atom);
        }
      }
      seed = old_seed;

      touched_indices[i] = 0; /* down-dating */
    }
  }

  result += countRingSizeFeatures(state, fp_counts, ncounts, nBonds, which_bits,
                                  exclude_atom);
  /* Set the color property to represent all different atom types */
  ap = state.atoms.data();
  for (int i = 0; i < nAtoms; i++, ap++) {
    const auto symbol = (*ap)->getSymbol();
    atomColor(state, *ap) =
        symbol == "*" ? 0
                      : PeriodicTable::getTable()->getAtomicNumber(symbol);
    if (atomColor(state, *ap) <= 1) {
      atomColor(state, *ap) = 0; /* ignore hydrogens */
    }
    /* mark special atom types */
    if (atomColor(state, *ap) > 115) {
      atomColor(state, *ap) = -1;
    }
    if ((*ap)->getSymbol() == "A") {
      atomColor(state, *ap) = -1;
    }
    if (i + 1 == exclude_atom) {
      atomColor(state, *ap) = 0;
    }
    if (atomColor(state, *ap) > 1 &&
        (!as_query || atomSubDescriptor(state, *ap) == SUB_AS_IS ||
         (atomSubDescriptor(state, *ap) != NONE &&
          atomSubDescriptor(state, *ap) != SUB_MORE &&
          atomSubDescriptor(state, *ap) == degree[i] + SUB_ONE - 1))) {
      atomColor(state, *ap) += 32 * degree[i];
    } else {
      atomColor(state, *ap) = 0;
    }
  }
  bp = state.bonds.data();
  for (int i = 0; i < nBonds; i++, bp++) {
    if (SINGLE <= bondType(state, *bp) && bondType(state, *bp) <= ANY_BOND) {
      bondColor(state, *bp) = 5;
    } else {
      bondColor(state, *bp) = 0;
    }
    if (bondEndpoint(state, (*bp)->getIdx(), 0) == exclude_atom) {
      bondColor(state, *bp) = 0;
    }
    if (bondEndpoint(state, (*bp)->getIdx(), 1) == exclude_atom) {
      bondColor(state, *bp) = 0;
    }
  }

  /*
   * add special bits for paths starting at atoms with >= 4 neighbours or
   * methyl atoms
   *
   * atom color represent degree bond color identical for all bonds to non-H
   */
  if (which_bits & USE_DEGREE_PATH) {
    seed = DEGREE_PATH_SEED;

    // set bits for degree paths starting with special carbon atoms
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      if (atomColor(state, *ap) <= 0) {
        continue;  // only process if degree defined
      }
      if (i + 1 == exclude_atom) {
        continue;
      }
      /* don't start on usual atoms */
      if (degree[i] <= 3 &&       // include if high degree
          atom_status[i] <= 2 &&  // include if ring fusion
          degree[i] != 1) {       // include if terminal atom
        continue;
      }
      /* only keep terminals if methyl */
      if (degree[i] == 1 && (*ap)->getSymbol() != "C") {
        continue;
      }
      // (*ap)->getSymbol(), i+1, degree[i], atom_status[i]);
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, atomColor(state, *ap));
      // (*ap)->getSymbol(), i+1, degree[i], atomSubDescriptor(state, *ap));
      result +=
          SetPathBitsRec(mol, state, nbp, fp_counts, ncounts, seed,
                         touched_indices, 1, 2, 4, /* path length 1 to 3 */
                         i, 0, -1,
                         // IGNORE_TERM_SYMBOL |
                         IGNORE_PATH_SYMBOL | PROCESS_CHAINS, exclude_atom);
      /* special CH fusion atoms */
      if (atom_status[i] > 2 && H_count[i + 1] >= 1) {
        seed = NEXT_SEED(seed, 219);
        result += SetPathBitsRec(
            mol, state, nbp, fp_counts, ncounts, seed, touched_indices, 1, 2, 5,
            i, 0, -1, IGNORE_PATH_SYMBOL | PROCESS_CHAINS, exclude_atom);
      }
      seed = old_seed;
      touched_indices[i] = 0; /* down-dating */
    }

    // set bits for degree paths starting with hetero atoms
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      if (i + 1 == exclude_atom) {
        continue;
      }
      const auto symbol = (*ap)->getSymbol();
      const auto tmp =
          symbol == "*" ? 0
                        : PeriodicTable::getTable()->getAtomicNumber(symbol);
      if (1 >= tmp || tmp >= 115) {
        continue;
      }
      if (tmp == 6) {
        continue;
      }
      if (tmp < 10) {
        continue;  // exclude common hetero atoms
      }
      touched_indices[i] = 1; /* updating */
      old_seed = seed;
      seed = NEXT_SEED(seed, tmp);
      // (*ap)->getSymbol(), i+1, degree[i], atomSubDescriptor(state, *ap));
      result += SetPathBitsRec(mol, state, nbp, fp_counts, ncounts, seed,
                               touched_indices, 1, 2, 2, i, 0, -1,
                               PROCESS_CHAINS, exclude_atom);
      if (0) {  // might overly populate complexes
        result += SetPathBitsRec(
            mol, state, nbp, fp_counts, ncounts, 101 + seed, touched_indices, 1,
            2, 2, i, 0, -1, FORCED_RING_PATH | PROCESS_CHAINS, exclude_atom);
      }
      seed = old_seed;
      touched_indices[i] = 0; /* down-dating */
    }
  }

  std::vector<std::vector<int>> length_matrix;
  if (which_bits & (USE_CLASS_SPIDERS | USE_FEATURE_PAIRS | USE_NON_SSS_BITS)) {
    /* Collect length_matrix */
    /* allocate storage length_matrix */
    length_matrix.assign(nAtoms, std::vector<int>(nAtoms, 0));
    for (int i = 0; i < nAtoms; i++) {
      touched_indices[i] = 0;
    }
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      if (i + 1 == exclude_atom) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      // if (false) fprintf(stderr, "starting path search at atom %d(%d)\n",
      // i+1, atomColor(state, *ap));
      SetPathLengthFlags(mol, state, touched_indices, i, 0, i,
                         12, /* path perception distance <= 12 */
                         length_matrix, nbp, exclude_atom);
      touched_indices[i] = 0; /* down-dating */
    }
  }

  /*
   * This screen class will catch non-linear fragments composed of rather
   * frequent linear sub-fragments
   */
  if (which_bits & (USE_CLASS_SPIDERS | USE_FEATURE_PAIRS)) {
    constexpr int MAX_SPIDER = 7;
    std::array<int, MAX_SPIDER + 1> csp3{};
    std::array<int, MAX_SPIDER + 1> hetero{};
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      const auto symbol = (*ap)->getSymbol();
      atomColor(state, *ap) =
          symbol == "*" ? 0
                        : PeriodicTable::getTable()->getAtomicNumber(symbol);
      if ((*ap)->getSymbol() == "H") {
        atomColor(state, *ap) = 0; /* ignore hydrogens */
      } else if ((*ap)->getSymbol() == "D") {
        atomColor(state, *ap) = 0; /* ignore hydrogens */
      } else if ((*ap)->getSymbol() == "T") {
        atomColor(state, *ap) = 0; /* ignore hydrogens */
      } else if ((*ap)->getSymbol() == "Q") {
        atomColor(state, *ap) = HETERO;
      } else if ((*ap)->getSymbol() == "A") {
        atomColor(state, *ap) = GENERIC;
      } else if ((*ap)->getSymbol() == "L") {
        atomColor(state, *ap) = GENERIC;
      } else if ((*ap)->getSymbol() == "C") {
        atomColor(state, *ap) =
            6; /* carbon second row elements are one class */
        if (cdegree[i] >= 3) {
          atomColor(state, *ap) = CSP3;
        }
      } else if (atomColor(state, *ap) > 1 && atomColor(state, *ap) < 115) {
        atomColor(state, *ap) = HETERO;
      } else { /* This could be R atoms or other odd things */
        atomColor(state, *ap) = 0;
      }
      if (i + 1 == exclude_atom) {
        atomColor(state, *ap) = 0;
      }
    }
    /*
     * Bond colors are already set OK, i.e. equal for A-H bonds
     */
    /* NOP */

    /* Now we start setting bits */
    ap = state.atoms.data();
    for (int i = 0; i < nAtoms; i++, ap++) {
      if (i + 1 == exclude_atom) {
        continue;
      }
      /* Spiders have at least three legs (;-) */
      if (degree[i] < 3) {
        continue;
      }
      /* Spider needs to be special atom or carbon */
      if (atomColor(state, *ap) != CSP3 && atomColor(state, *ap) != 6) {
        continue;
      }
      touched_indices[i] = 1; /* updating */
      std::fill(csp3.begin(), csp3.end(), 0);
      std::fill(hetero.begin(), hetero.end(), 0);
      if (which_bits & USE_CLASS_SPIDERS) {
        SpecialNeighboursRec(mol, state, touched_indices, 1, i, MAX_SPIDER,
                             csp3.data(), hetero.data(), nbp, exclude_atom);
      }
      touched_indices[i] = 0; /* down-dating */

      /* set bits for spiders with one CSP3 atom and two heteros */
      if (which_bits & USE_CLASS_SPIDERS) {
        for (int j = 1; j <= MAX_SPIDER; j++) {
          if (csp3[j] == 0) {
            continue;
          }
          seed = CLASS_SPIDER_SEED;
          if (atomColor(state, *ap) == HETERO) {
            // seed = CLASS_SPIDER_SEED+HETERO*8+CSP3*11;
            seed = NEXT_SEED(seed, HETERO * 8);
            seed = NEXT_SEED(seed, CSP3 * 11);
          } else {
            // seed = CLASS_SPIDER_SEED+6*8+CSP3*11;
            seed = NEXT_SEED(seed, 6 * 8);
            seed = NEXT_SEED(seed, CSP3 * 11);
          }
          for (int j1 = 1; j1 <= MAX_SPIDER; j1++) {
            const auto tmp1 = hetero[j1];
            if (tmp1 <= 0) {
              continue;
            }
            for (int j2 = j1; j2 <= MAX_SPIDER; j2++) {
              auto tmp2 = hetero[j2];
              if (j2 == j1) {
                tmp2--; /* consumed in outer loop */
              }
              if (tmp2 <= 0) {
                continue;
              }
              old_seed = seed;
              seed = NEXT_SEED(seed, j);
              seed = NEXT_SEED(seed, j1 + j2);
              ADD_BIT(fp_counts, ncounts, seed);
              result++;
              seed = old_seed;
            }
          }
        }
      }

      /* don't hetero-spider normal carbons */
      if (atomColor(state, *ap) != CSP3) {
        continue;
      }
      /* set bits for spiders with three defined HETERO atoms */
      if (which_bits & USE_CLASS_SPIDERS) {
        for (int j = 1; j <= MAX_SPIDER; j++) {
          if (hetero[j] == 0) {
            continue;
          }
          seed = CLASS_SPIDER_SEED;
          if (atomColor(state, *ap) == HETERO) {
            // seed = CLASS_SPIDER_SEED+HETERO*8+HETERO*11;
            seed = NEXT_SEED(seed, HETERO * 8);
            seed = NEXT_SEED(seed, HETERO * 11);
          } else {
            // seed = CLASS_SPIDER_SEED+6*8+HETERO*11;
            seed = NEXT_SEED(seed, 6 * 8);
            seed = NEXT_SEED(seed, HETERO * 11);
          }
          for (int j1 = j; j1 <= MAX_SPIDER; j1++) {
            auto tmp1 = hetero[j1];
            if (j1 == j) {
              tmp1--; /* we've consumed this one in outer loop */
            }
            if (tmp1 <= 0) {
              continue;
            }
            for (int j2 = j1; j2 <= MAX_SPIDER; j2++) {
              auto tmp2 = hetero[j2];
              if (j2 == j) {
                tmp2--; /* consumed in outer loop */
              }
              if (j2 == j1) {
                tmp2--; /* consumed in outer loop */
              }
              if (tmp2 <= 0) {
                continue;
              }
              old_seed = seed;
              seed = NEXT_SEED(seed, j * j1 * j2);
              ADD_BIT(fp_counts, ncounts, seed);
              result++;
              // Additional bits for quarternary centers
              if (degree[i] > 3) {
                ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, 501));
                result++;
              }
              seed = old_seed;
            }
          }
        }
      }
    }

    /**
     * Collect bits that represent feature/path_length/feature triples.
     */
    result += countFeaturePairFeatures(mol, state, fp_counts, ncounts, nAtoms,
                                       atom_status, degree, cdegree,
                                       length_matrix, which_bits, exclude_atom);
  }
  result += countScaffoldFeatures(state, nbp, fp_counts, ncounts, atom_status,
                                  bond_status, degree, length_matrix, nAtoms,
                                  which_bits, as_query, exclude_atom);
  if (which_bits & (USE_CLASS_SPIDERS | USE_FEATURE_PAIRS | USE_NON_SSS_BITS)) {
    /* de-allocate length_matrix */
  }

  // seed = 0;
  // for (i=0; i<nbytes; i++)
  //    seed = (0xFF&fingerprint[i]) ^ ((seed>>8) | ((seed&0xFF)<<16));
  //         CountBits(fingerprint, nbytes), seed, result);
  // FORTIFY    fortifyTest(total_bytes_allocated,
  // "SetFingerprintCountsWithFocus");
  return (result);
}

// Enumerates the rings through ring atoms/bonds (up to max_size) and sets the
// ring size flags (bit k: member of a ring of size k; bit 0: any ring)
void markRingsRecursive(const ROMol &mol, AvalonState &state,
                        std::vector<int> &touchedAtoms,
                        std::vector<int> &touchedBonds, int startIndex,
                        int pathLength, int currentIndex, int maxSize,
                        const std::vector<neighbourhood_t> &nbp) {
  const auto nAtoms = static_cast<int>(mol.getNumAtoms());
  const auto nBonds = static_cast<int>(mol.getNumBonds());
  for (int i = 0; i < nbp[currentIndex].n_ligands; i++) {
    int ai = nbp[currentIndex].atoms[i];
    if (ai < startIndex) {
      continue;
    }
    if (ai == startIndex) {
      if (pathLength < 3) {
        continue;
      }
      for (int j = 0; j < nAtoms; j++) {
        if (touchedAtoms[j]) {
          state.atomRingFlags[j] |= (1 << pathLength);
        }
      }
      for (int j = 0; j < nBonds; j++) {
        if (touchedBonds[j]) {
          state.bondRingFlags[j] |= (1 << pathLength);
        }
      }
      continue;
    }
    if (touchedAtoms[ai]) {
      continue;
    }
    if (pathLength + 1 > maxSize) {
      continue;
    }
    if (state.atomRingFlags[ai] == 0) {
      continue;
    }
    int bi = nbp[currentIndex].bonds[i];
    if (state.bondRingFlags[bi] == 0) {
      continue;
    }
    touchedAtoms[ai] = 1;
    touchedBonds[bi] = 1;
    markRingsRecursive(mol, state, touchedAtoms, touchedBonds, startIndex,
                       pathLength + 1, ai, maxSize, nbp);
    touchedAtoms[ai] = 0;
    touchedBonds[bi] = 0;
  }
}

}  // namespace

namespace {
using BondSet = std::vector<char>;

const int kSingle = 1;
const int kDouble = 2;
const int kTriple = 3;
const int kAromatic = 4;

// The ring sets used by the perception code: the (symmetrized) SSSR rings
// which only use usable bonds, followed by the rings obtainable by XORing two
// base rings which have exactly one connected path in common.
struct RingSets {
  std::vector<BondSet> base;
  std::vector<BondSet> pairs;
};

int cardinality(const BondSet &bs) {
  int res = 0;
  for (auto c : bs) {
    res += c != 0;
  }
  return res;
}

RingSets findRingSets(const ROMol &mol, const std::vector<char> &usableBond) {
  RingSets res;
  const auto nBonds = mol.getNumBonds();
  for (const auto &ring : mol.getRingInfo()->bondRings()) {
    bool ok = true;
    for (auto bi : ring) {
      ok = ok && usableBond[bi];
    }
    if (!ok) {
      continue;
    }
    BondSet bs(nBonds, 0);
    for (auto bi : ring) {
      bs[bi] = 1;
    }
    res.base.push_back(std::move(bs));
  }
  std::vector<int> touched(mol.getNumAtoms());
  for (size_t i = 0; i < res.base.size(); ++i) {
    for (size_t j = i + 1; j < res.base.size(); ++j) {
      std::fill(touched.begin(), touched.end(), 0);
      bool overlap = false;
      for (unsigned int b = 0; b < nBonds; ++b) {
        if (res.base[i][b] && res.base[j][b]) {
          overlap = true;
          const auto bond = mol.getBondWithIdx(b);
          ++touched[bond->getBeginAtomIdx()];
          ++touched[bond->getEndAtomIdx()];
        }
      }
      if (!overlap) {
        continue;
      }
      if (std::count(touched.begin(), touched.end(), 1) == 2) {
        BondSet bs(nBonds, 0);
        for (unsigned int b = 0; b < nBonds; ++b) {
          bs[b] = res.base[i][b] != res.base[j][b];
        }
        res.pairs.push_back(std::move(bs));
      }
    }
  }
  return res;
}

// port of PerceiveAromaticBonds(): six-rings (and fused rings with 4n+2
// bonds) which only consist of sp2 atoms are aromatic
void perceiveAromaticBonds(const ROMol &mol, AvalonState &state) {
  const int nBonds = mol.getNumBonds();
  const int nAtoms = mol.getNumAtoms();
  std::vector<char> allBonds(nBonds, 1);
  auto rs = findRingSets(mol, allBonds);
  std::vector<char> inRing(nBonds, 0);
  for (const auto &r : rs.base) {
    for (int b = 0; b < nBonds; ++b) {
      inRing[b] |= r[b];
    }
  }
  std::vector<const BondSet *> rings;
  for (const auto &r : rs.pairs) {
    rings.push_back(&r);
  }
  for (const auto &r : rs.base) {
    rings.push_back(&r);
  }
  std::vector<int> spCount(nAtoms + 1);
  bool changed;
  do {
    changed = false;
    for (const auto ring : rings) {
      std::fill(spCount.begin(), spCount.end(), 0);
      int ndouble = 0, nsingle = 0;
      for (int b = 0; b < nBonds; ++b) {
        if ((*ring)[b] && state.bondTypes[b] == kSingle) {
          ++nsingle;
        }
      }
      bool isCumulene = false;
      for (int b = 0; b < nBonds; ++b) {
        if (!(*ring)[b]) {
          continue;
        }
        if (state.bondTypes[b] == kDouble) {
          ++ndouble;
          for (int k = 0; k < 2; ++k) {
            if (++spCount[bondEndpoint(state, b, k)] > 1) {
              isCumulene = true;
            }
          }
        } else if (state.bondTypes[b] == kTriple) {
          isCumulene = true;
        }
      }
      for (int b = 0; b < nBonds; ++b) {
        if ((*ring)[b] && state.bondTypes[b] == kAromatic) {
          for (int k = 0; k < 2; ++k) {
            if (spCount[bondEndpoint(state, b, k)] == 0) {
              ++spCount[bondEndpoint(state, b, k)];
            }
          }
        }
      }
      bool isAromatic = !isCumulene && ((cardinality(*ring) - 2) % 4) == 0;
      for (int b = 0; b < nBonds; ++b) {
        if ((*ring)[b] && (spCount[bondEndpoint(state, b, 0)] != 1 ||
                           spCount[bondEndpoint(state, b, 1)] != 1)) {
          isAromatic = false;
        }
      }
      if (isAromatic && (ndouble > 0 || nsingle > 0)) {
        for (int b = 0; b < nBonds; ++b) {
          if ((*ring)[b] && inRing[b] && state.bondTypes[b] != kAromatic) {
            changed = true;
            state.bondTypes[b] = kAromatic;
          }
        }
      }
    }
  } while (changed);
}

// port of PerceiveDYAromaticity(): Daylight-like aromaticity
void perceiveDYAromaticity(const ROMol &mol, AvalonState &state,
                           const std::vector<neighbourhood_t> &nbp) {
  const int nAtoms = mol.getNumAtoms();
  const int nBonds = mol.getNumBonds();
  std::vector<char> allBonds(nBonds, 1);
  auto allRings = findRingSets(mol, allBonds);
  if (allRings.base.empty()) {
    return;
  }
  std::vector<char> bondInRing(nBonds, 0);
  for (const auto &r : allRings.base) {
    for (int b = 0; b < nBonds; ++b) {
      bondInRing[b] |= r[b];
    }
  }
  std::vector<char> atomInRing(nAtoms, 0);
  std::vector<char> candidate(nAtoms, 0);
  for (int i = 0; i < nAtoms; ++i) {
    candidate[i] = mol.getAtomWithIdx(i)->getSymbol() != "C";
  }
  for (int b = 0; b < nBonds; ++b) {
    if (state.bondTypes[b] > kSingle && state.bondTypes[b] != kTriple) {
      candidate[bondEndpoint(state, b, 0) - 1] = 1;
      candidate[bondEndpoint(state, b, 1) - 1] = 1;
    }
  }
  std::vector<char> usable(nBonds, 0);
  for (int b = 0; b < nBonds; ++b) {
    if (!bondInRing[b]) {
      continue;
    }
    atomInRing[bondEndpoint(state, b, 0) - 1] = 1;
    atomInRing[bondEndpoint(state, b, 1) - 1] = 1;
    usable[b] = candidate[bondEndpoint(state, b, 0) - 1] &&
                candidate[bondEndpoint(state, b, 1) - 1];
  }
  auto rs = findRingSets(mol, usable);
  if (rs.base.empty()) {
    return;
  }
  std::vector<BondSet> candidates;
  for (const auto &r : rs.base) {
    candidates.push_back(r);
  }
  for (const auto &r : rs.pairs) {
    candidates.push_back(r);
  }
  // fused pairs
  const size_t nRings = candidates.size();
  for (size_t i = 0; i < nRings; ++i) {
    for (size_t j = i + 1; j < nRings; ++j) {
      bool overlap = false;
      BondSet x(nBonds, 0);
      for (int b = 0; b < nBonds; ++b) {
        overlap = overlap || (candidates[i][b] && candidates[j][b]);
        x[b] = candidates[i][b] != candidates[j][b];
      }
      if (overlap && cardinality(x) == cardinality(candidates[i]) +
                                           cardinality(candidates[j]) - 2) {
        candidates.push_back(std::move(x));
      }
    }
  }

  const auto isSym = [&](int ai, const std::vector<const char *> &list) {
    return AtomSymbolMatch(mol.getAtomWithIdx(ai)->getSymbol(), list);
  };
  bool changed;
  do {
    changed = false;
    for (const auto &ring : candidates) {
      int npi = 0;
      bool conjugated = true;
      for (int i = 0; i < nAtoms; ++i) {
        if (!atomInRing[i]) {
          continue;
        }
        int inRingDouble = 0, inRingAromatic = 0;
        bool exoPull = false, isInRing = false;
        int localPi = 0;
        for (int j = 0; j < nbp[i].n_ligands; ++j) {
          const int bi = nbp[i].bonds[j];
          const int ai = nbp[i].atoms[j];
          if (bondInRing[bi] && ring[bi]) {
            isInRing = true;
            if (state.bondTypes[bi] == kAromatic) {
              ++inRingAromatic;
            }
            if (state.bondTypes[bi] == kDouble) {
              ++inRingDouble;
            }
          } else if (state.bondTypes[bi] == kDouble) {
            if (mol.getAtomWithIdx(i)->getSymbol() == "C" && !bondInRing[bi] &&
                isSym(ai, {"O", "S", "P", "N", "L"})) {
              exoPull = true;
            }
          }
        }
        if (!isInRing) {
          continue;
        }
        if ((inRingAromatic >= 1 || inRingDouble == 1) &&
            (isSym(i, {"C", "N", "A", "*"}) ||
             mol.getAtomWithIdx(i)->getSymbol() == "L")) {
          localPi = 1;
        } else if (inRingAromatic == 0 && inRingDouble == 0 &&
                   mol.getAtomWithIdx(i)->getFormalCharge() == 0 &&
                   isSym(i, {"N", "S", "O"})) {
          localPi = 2;
        } else if (inRingAromatic == 0 && inRingDouble == 0 &&
                   mol.getAtomWithIdx(i)->getFormalCharge() == 0 && exoPull &&
                   mol.getAtomWithIdx(i)->getSymbol() == "C") {
          localPi = 0;
        } else {
          conjugated = false;
        }
        if (mol.getAtomWithIdx(i)->getFormalCharge() < 0 &&
            isSym(i, {"C", "N"})) {
          conjugated = false;
        }
        npi += localPi;
      }
      if (!conjugated || npi % 4 != 2) {
        continue;
      }
      for (int b = 0; b < nBonds; ++b) {
        if (bondInRing[b] && ring[b] && state.bondTypes[b] != kAromatic) {
          state.bondTypes[b] = kAromatic;
          changed = true;
        }
      }
    }
  } while (changed);

  std::fill(candidate.begin(), candidate.end(), 0);
  for (int b = 0; b < nBonds; ++b) {
    if (state.bondTypes[b] == kAromatic) {
      candidate[bondEndpoint(state, b, 0) - 1] =
          candidate[bondEndpoint(state, b, 1) - 1] = 1;
    }
  }
  for (int b = 0; b < nBonds; ++b) {
    if (candidate[bondEndpoint(state, b, 0) - 1] &&
        candidate[bondEndpoint(state, b, 1) - 1] && bondInRing[b] &&
        state.bondTypes[b] == kSingle) {
      state.bondTypes[b] = kAromatic;
    }
  }
}

// port of the carbon part of GuessHCountsFromSubstitution() as it applies to
// molecules without substitution-count queries
void guessSubstitution(const ROMol &mol, AvalonState &state,
                       const std::vector<neighbourhood_t> &nbp) {
  for (size_t i = 0; i < mol.getNumAtoms(); ++i) {
    if (mol.getAtomWithIdx(i)->getFormalCharge() != 0 ||
        mol.getAtomWithIdx(i)->getNumRadicalElectrons() != 0 ||
        mol.getAtomWithIdx(i)->getSymbol() != "C") {
      continue;
    }
    int nsingle = 0, ndouble = 0, ntriple = 0, naromatic = 0, nother = 0;
    for (int j = 0; j < nbp[i].n_ligands; ++j) {
      switch (state.bondTypes[nbp[i].bonds[j]]) {
        case kSingle:
          ++nsingle;
          break;
        case kDouble:
          ++ndouble;
          break;
        case kAromatic:
          ++naromatic;
          break;
        case kTriple:
          ++ntriple;
          break;
        default:
          ++nother;
      }
    }
    if (nother > 0) {
      continue;
    }
    if (naromatic == 2) {
      ++nsingle;
      ++ndouble;
      naromatic = 0;
    }
    if (naromatic > 0) {
      continue;
    }
    if (nsingle + 2 * ndouble + 3 * ntriple == 4) {
      state.atomSubDescriptors[i] = SUB_AS_IS;
    }
  }
}
}  // namespace

std::vector<std::uint32_t> getAvalonCounts(const ROMol &mol,
                                           std::uint32_t fpSize,
                                           std::uint32_t bitFlags, bool isQuery,
                                           int focusAtom,
                                           bool accumulateAsQuery) {
  PRECONDITION(fpSize > 0, "fpSize must be > 0");
  PRECONDITION(
      focusAtom < 0 || static_cast<unsigned>(focusAtom) < mol.getNumAtoms(),
      "focusAtom out of range");
  std::vector<std::uint32_t> res(fpSize, 0);
  const int nAtoms = mol.getNumAtoms();
  if (!nAtoms) {
    return res;
  }

  // Avalon does its own aromaticity perception, so work with a kekulized copy
  auto tmol = std::make_unique<RWMol>(mol);
  try {
    MolOps::Kekulize(*tmol, true);
  } catch (const MolSanitizeException &) {
    // keep the aromatic bonds as they are
  }
  if (!tmol->getRingInfo()->isInitialized()) {
    MolOps::fastFindRings(*tmol);
  }
  if (!isQuery) {
    for (const auto atom : tmol->atoms()) {
      if (atom->needsUpdatePropertyCache()) {
        tmol->updatePropertyCache(false);
        break;
      }
    }
  }
  const auto *lmol = tmol.get();
  const auto &ri = *lmol->getRingInfo();
  const int nBonds = lmol->getNumBonds();

  AvalonState state(*lmol);
  std::vector<neighbourhood_t> nbp(nAtoms);
  std::vector<int> explicitH(nAtoms + 1, 0);
  std::vector<int> atomStatus(nAtoms, 0);
  std::vector<int> bondStatus(nBonds, 0);

  for (const auto bond : lmol->bonds()) {
    const int bi = bond->getIdx();
    const int bai = bond->getBeginAtomIdx();
    const int eai = bond->getEndAtomIdx();
    nbp[bai].atoms.push_back(eai);
    nbp[bai].bonds.push_back(bi);
    nbp[bai].n_ligands++;
    nbp[eai].atoms.push_back(bai);
    nbp[eai].bonds.push_back(bi);
    nbp[eai].n_ligands++;
    // explicit hydrogens in the graph count towards the H count
    if (lmol->getAtomWithIdx(eai)->getAtomicNum() == 1) {
      explicitH[bai + 1]++;
    } else if (lmol->getAtomWithIdx(bai)->getAtomicNum() == 1) {
      explicitH[eai + 1]++;
    }
    bondStatus[bi] = ri.numBondRings(bi);
    if (bondStatus[bi]) {
      atomStatus[bai]++;
      atomStatus[eai]++;
    }
  }
  const auto kekuleBondTypes = state.bondTypes;

  for (int i = 0; i < nAtoms; ++i) {
    state.atomRingFlags[i] = atomStatus[i] > 0 ? 1 : 0;
  }
  for (int i = 0; i < nBonds; ++i) {
    state.bondRingFlags[i] = bondStatus[i] > 0 ? 1 : 0;
  }
  std::vector<int> touchedAtoms(nAtoms, 0);
  std::vector<int> touchedBonds(nBonds, 0);
  for (int i = 0; i < nAtoms; ++i) {
    if (state.atomRingFlags[i] == 0) {
      continue;
    }
    touchedAtoms[i] = 1;
    markRingsRecursive(*lmol, state, touchedAtoms, touchedBonds, i, 1, i, 14,
                       nbp);
    touchedAtoms[i] = 0;
  }

  // The Avalon wrapper in AvalonTools.cpp runs the algorithm twice for
  // non-query molecules: once with Avalon's simple aromaticity model and once
  // with Daylight-like aromaticity, accumulating the results.
  std::vector<int> counts(fpSize, 0);
  const int nPasses = isQuery ? 1 : 2;
  for (int pass = 0; pass < nPasses; ++pass) {
    const bool dy = pass == 1;
    const bool queryMode = dy ? accumulateAsQuery : isQuery;
    for (int i = 0; i < nBonds; ++i) {
      state.bondTypes[i] = kekuleBondTypes[i];
    }
    auto hCount = explicitH;
    std::fill(state.atomSubDescriptors.begin(), state.atomSubDescriptors.end(),
              0);
    if (!queryMode) {
      for (const auto atom : lmol->atoms()) {
        hCount[atom->getIdx() + 1] += atom->getTotalNumHs();
      }
    } else {
      guessSubstitution(*lmol, state, nbp);
    }
    if (dy) {
      perceiveDYAromaticity(*lmol, state, nbp);
    } else {
      perceiveAromaticBonds(*lmol, state);
    }
    CountFingerprintPatterns(
        *lmol, state, nbp, hCount.data(), atomStatus.data(), bondStatus.data(),
        counts.data(), fpSize, bitFlags, queryMode, focusAtom + 1);
  }
  for (unsigned int i = 0; i < fpSize; ++i) {
    res[i] = counts[i] > 0 ? counts[i] : 0;
  }
  return res;
}

// ---- FingerprintGenerator integration ----

AvalonArguments::AvalonArguments(std::uint32_t bitFlags, bool isQuery,
                                 bool countSimulation,
                                 const std::vector<std::uint32_t> countBounds,
                                 std::uint32_t fpSize, bool accumulateAsQuery)
    : FingerprintArguments(countSimulation, countBounds, fpSize, 1, false),
      d_bitFlags(bitFlags),
      df_isQuery(isQuery),
      df_accumulateAsQuery(accumulateAsQuery) {
  PRECONDITION(fpSize > 0, "fpSize must be > 0");
}

std::string AvalonArguments::infoString() const {
  return "AvalonArguments bitFlags=" + std::to_string(d_bitFlags) +
         " isQuery=" + std::to_string(df_isQuery) +
         " accumulateAsQuery=" + std::to_string(df_accumulateAsQuery) + " " +
         commonArgumentsString();
}

void AvalonArguments::toJSON(boost::property_tree::ptree &pt) const {
  pt.put("type", "AvalonArguments");
  pt.put("bitFlags", d_bitFlags);
  pt.put("isQuery", df_isQuery);
  pt.put("accumulateAsQuery", df_accumulateAsQuery);
  FingerprintArguments::toJSON(pt);
}

void AvalonArguments::fromJSON(const boost::property_tree::ptree &pt) {
  d_bitFlags = pt.get<std::uint32_t>("bitFlags", d_bitFlags);
  df_isQuery = pt.get<bool>("isQuery", df_isQuery);
  df_accumulateAsQuery =
      pt.get<bool>("accumulateAsQuery", df_accumulateAsQuery);
  FingerprintArguments::fromJSON(pt);
}

template <typename OutputType>
OutputType AvalonAtomEnv<OutputType>::getBitId(
    FingerprintArguments *, const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *, AdditionalOutput *, const bool,
    const std::uint64_t) const {
  return d_bitId;
}

template <typename OutputType>
void AvalonAtomEnv<OutputType>::updateAdditionalOutput(AdditionalOutput *,
                                                       std::uint64_t) const {}

template <typename OutputType>
std::vector<AtomEnvironment<OutputType> *>
AvalonEnvGenerator<OutputType>::getEnvironments(
    const ROMol &mol, FingerprintArguments *arguments,
    const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *ignoreAtoms, int,
    const AdditionalOutput *, const std::vector<std::uint32_t> *,
    const std::vector<std::uint32_t> *, bool) const {
  PRECONDITION(arguments, "bad arguments");
  if ((fromAtoms && !fromAtoms->empty()) ||
      (ignoreAtoms && !ignoreAtoms->empty())) {
    throw ValueErrorException(
        "fromAtoms and ignoreAtoms are not supported by the Avalon "
        "fingerprint generator");
  }
  auto *args = dynamic_cast<AvalonArguments *>(arguments);
  PRECONDITION(args, "arguments must be AvalonArguments");

  // the bit ids are the slots of the (hashed) Avalon count vector; each slot
  // contributes one environment per count
  auto counts =
      getAvalonCounts(mol, args->d_fpSize, args->d_bitFlags, args->df_isQuery,
                      -1, args->df_accumulateAsQuery);
  std::vector<AtomEnvironment<OutputType> *> res;
  for (std::uint32_t i = 0; i < counts.size(); ++i) {
    for (std::uint32_t j = 0; j < counts[i]; ++j) {
      res.push_back(new AvalonAtomEnv<OutputType>(i));
    }
  }
  return res;
}

template <typename OutputType>
std::string AvalonEnvGenerator<OutputType>::infoString() const {
  return "AvalonEnvGenerator fpSize=" + std::to_string(d_fpSize);
}

template <typename OutputType>
void AvalonEnvGenerator<OutputType>::toJSON(
    boost::property_tree::ptree &pt) const {
  pt.put("type", "AvalonEnvGenerator");
  pt.put("fpSize", d_fpSize);
  AtomEnvironmentGenerator<OutputType>::toJSON(pt);
}

template <typename OutputType>
void AvalonEnvGenerator<OutputType>::fromJSON(
    const boost::property_tree::ptree &pt) {
  d_fpSize = pt.get<std::uint32_t>("fpSize", d_fpSize);
  AtomEnvironmentGenerator<OutputType>::fromJSON(pt);
}

template <typename OutputType>
OutputType AvalonEnvGenerator<OutputType>::getResultSize() const {
  return d_fpSize;
}

template <typename OutputType>
FingerprintGenerator<OutputType> *getAvalonGenerator(
    const AvalonArguments &args) {
  return new FingerprintGenerator<OutputType>(
      new AvalonEnvGenerator<OutputType>(args.d_fpSize),
      new AvalonArguments(args), nullptr, nullptr, false, false);
}

template <typename OutputType>
FingerprintGenerator<OutputType> *getAvalonGenerator(std::uint32_t fpSize,
                                                     std::uint32_t bitFlags,
                                                     bool isQuery,
                                                     bool accumulateAsQuery) {
  AvalonArguments args(bitFlags, isQuery, false, {1, 2, 4, 8}, fpSize,
                       accumulateAsQuery);
  return getAvalonGenerator<OutputType>(args);
}

template class AvalonAtomEnv<std::uint32_t>;
template class AvalonAtomEnv<std::uint64_t>;
template class AvalonEnvGenerator<std::uint32_t>;
template class AvalonEnvGenerator<std::uint64_t>;

template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint32_t> *
getAvalonGenerator(const AvalonArguments &);
template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint64_t> *
getAvalonGenerator(const AvalonArguments &);
template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint32_t> *
getAvalonGenerator(std::uint32_t, std::uint32_t, bool, bool);
template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint64_t> *
getAvalonGenerator(std::uint32_t, std::uint32_t, bool, bool);

}  // namespace AvalonFP
}  // namespace RDKit
