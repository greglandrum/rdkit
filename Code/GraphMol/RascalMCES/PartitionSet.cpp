//
// Copyright (C) David Cosgrove 2023
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>

#include "PartitionSet.h"

namespace RDKit {

namespace RascalMCES {
PartitionSet::PartitionSet(const std::vector<std::vector<char>> &modProd,
                           const std::vector<std::pair<int, int>> &vtxPairs,
                           const std::vector<unsigned int> &vtx1Labels,
                           const std::vector<unsigned int> &vtx2Labels,
                           unsigned int lowerBound)
    : d_ModProd(new std::vector<std::vector<char>>(modProd)),
      d_VtxPairs(new std::vector<std::pair<int, int>>(vtxPairs)),
      d_vtx1Labels(new std::vector<unsigned int>(vtx1Labels)),
      d_vtx2Labels(new std::vector<unsigned int>(vtx2Labels)) {
  d_vtx1Counts = std::vector<int>(d_vtx1Labels->size(), 0);
  d_vtx2Counts = std::vector<int>(d_vtx2Labels->size(), 0);
  int firstVtx = -1;
  // Clearly, a vertex in one of the line graphs can only match one vertex
  // in the other.  Thus, the initial partitions can be set up so that
  // all vertices in a partition have the same vertex in the first
  // line graph.
  for (size_t i = 0; i < vtxPairs.size(); ++i) {
    auto &vp = vtxPairs[i];
    if (vp.first != firstVtx) {
      d_parts.push_back(std::vector<unsigned int>());
      d_parts.back().push_back(i);
      firstVtx = vp.first;
    } else {
      d_parts.back().push_back(i);
    }
    d_vtx1Counts[vp.first]++;
    d_vtx2Counts[vp.second]++;
  }
  if (d_parts.empty()) {
    return;
  }
  // Now sort the partitions by size.  This means that the vertices at the
  // top of the partition set, above the lowerBound (or Pex as Raymond
  // calls it in the paper), are the ones that match the least number of
  // vertices in the other line graph.  This has a dramatic effect on the
  // speed compared with other things tried.  I think it is what Raymond
  // means when he says "Perform an initial partitioning of the vertices...
  // using the labeled edge projection procedure."
  sort_partitions();
  // Now reassign vertices from above Pex to below it if possible.
  // This also improves the speed of finding a large clique early.
  // A vertex is moved to a partition where it isn't connected to a vertex
  // in the modular product graph that is in the partition.
  for (size_t i = d_parts.size() - 1; i > lowerBound; --i) {
    bool reassigned = false;
    for (auto &iv : d_parts[i]) {
      for (size_t k = 0; k <= lowerBound; ++k) {
        bool conn = false;
        for (auto kv : d_parts[k]) {
          if (modProd[iv][kv]) {
            conn = true;
            break;
          }
        }
        if (!conn) {
          d_parts[k].push_back(iv);
          iv = std::numeric_limits<unsigned int>::max();
          reassigned = true;
          break;
        }
      }
    }
    if (reassigned) {
      d_parts[i].erase(std::remove(d_parts[i].begin(), d_parts[i].end(),
                                   std::numeric_limits<unsigned int>::max()),
                       d_parts[i].end());
    }
  }
  d_parts.erase(std::remove_if(d_parts.begin(), d_parts.end(),
                               [](const std::vector<unsigned int> &v) {
                                 return v.empty();
                               }),
                d_parts.end());
  // Sort again, to make sure the large partitions are dealt with as late as
  // possible.
  sort_partitions();

  // Get the info together for the upper bound calculation.
  calc_vtx_type_counts();
}

int PartitionSet::upper_bound() {
  int upperBound = 0;
  for (size_t i = 0; i < d_vtx1TypeCounts.size(); ++i) {
    upperBound += std::min(d_vtx1TypeCounts[i], d_vtx2TypeCounts[i]);
  }
  return upperBound;
}

void PartitionSet::print_partitions(std::ostream &os) const {
  for (size_t i = 0; i < d_parts.size(); ++i) {
    os << i << " :: " << d_parts[i].size() << " ::";
    for (auto &mem : d_parts[i]) {
      os << " " << mem << " (" << (*d_VtxPairs)[mem].first << ","
         << (*d_VtxPairs)[mem].second << ")";
    }
    os << std::endl;
  }
  os << "vtx1_counts :";
  for (auto vc : d_vtx1Counts) {
    os << " " << vc;
  }
  os << std::endl;
  os << "vtx2_counts :";
  for (auto vc : d_vtx2Counts) {
    os << " " << vc;
  }
  os << std::endl;
}

unsigned int PartitionSet::last_vertex() const {
  if (d_parts.empty()) {
    return std::numeric_limits<unsigned int>::max();
  }
  return d_parts.back().back();
}

unsigned int PartitionSet::pop_last_vertex() {
  if (d_parts.empty()) {
    throw std::runtime_error("PartitionSet set is empty.");
  }
  unsigned int ret_val = d_parts.back().back();
  d_parts.back().pop_back();
  if (d_parts.back().empty()) {
    d_parts.pop_back();
  }
  decrement_vertex_counts(ret_val);
  return ret_val;
}

void PartitionSet::prune_vertices(unsigned int vtx_num) {
  for (auto &part : d_parts) {
    size_t i = 0;
    while (i < part.size()) {
      if (!(*d_ModProd)[part[i]][vtx_num]) {
        decrement_vertex_counts(part[i]);
        part[i] = part.back();
        part.pop_back();
      } else {
        ++i;
      }
    }
  }
  d_parts.erase(std::remove_if(d_parts.begin(), d_parts.end(),
                               [](const std::vector<unsigned int> &v) {
                                 return v.empty();
                               }),
                d_parts.end());
  sort_partitions();
}

void PartitionSet::add_vertex(unsigned int vtxNum) {
  // add vtxNum to the first partition where it isn't connected in
  // d_ModProd to any of the existing members.  Puts it in a new
  // partition if necessary.
  d_vtx1Counts[(*d_VtxPairs)[vtxNum].first]++;
  d_vtx2Counts[(*d_VtxPairs)[vtxNum].second]++;
  for (auto &part : d_parts) {
    bool conn = false;
    for (auto &mem : part) {
      if ((*d_ModProd)[mem][vtxNum]) {
        conn = true;
        break;
      }
    }
    if (!conn) {
      part.push_back(vtxNum);
      return;
    }
  }
  d_parts.push_back(std::vector<unsigned int>(1, vtxNum));
}

void PartitionSet::sort_partitions() {
  std::sort(d_parts.begin(), d_parts.end(),
            [](const std::vector<unsigned int> &v1,
               const std::vector<unsigned int> &v2) {
              return v1.size() > v2.size();
            });
}

void PartitionSet::calc_vtx_type_counts() {
  auto doIt = [](unsigned int maxLabel, const std::vector<int> &vtxCounts,
                 const std::vector<unsigned int> &vtxLabels,
                 std::vector<int> &vtxTypeCounts) -> void {
    vtxTypeCounts = std::vector<int>(maxLabel + 1, 0);
    for (size_t i = 0; i < vtxCounts.size(); ++i) {
      if (vtxCounts[i]) {
        ++vtxTypeCounts[vtxLabels[i]];
      }
    }
  };

  unsigned int max_label = 0;
  max_label =
      std::max(*std::max_element(d_vtx1Labels->begin(), d_vtx1Labels->end()),
               *std::max_element(d_vtx2Labels->begin(), d_vtx2Labels->end()));
  doIt(max_label, d_vtx1Counts, *d_vtx1Labels, d_vtx1TypeCounts);
  doIt(max_label, d_vtx2Counts, *d_vtx2Labels, d_vtx2TypeCounts);
}

void PartitionSet::decrement_vertex_counts(int vtxNum) {
  --d_vtx1Counts[(*d_VtxPairs)[vtxNum].first];
  if (!d_vtx1Counts[(*d_VtxPairs)[vtxNum].first]) {
    --d_vtx1TypeCounts[(*d_vtx1Labels)[(*d_VtxPairs)[vtxNum].first]];
  }
  --d_vtx2Counts[(*d_VtxPairs)[vtxNum].second];
  if (!d_vtx2Counts[(*d_VtxPairs)[vtxNum].second]) {
    --d_vtx2TypeCounts[(*d_vtx2Labels)[(*d_VtxPairs)[vtxNum].second]];
  }
}

}  // namespace RascalMCES
}  // namespace RDKit
