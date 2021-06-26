#include "permutation.hpp"
#include <algorithm>
#include <cassert>

Permutation::Permutation(const Permutation &pi) {
  genes = unique_ptr<vector<Gene>>(new vector<Gene>(*pi.genes));
  intergenic_regions =
      unique_ptr<vector<IR>>(new vector<IR>(*pi.intergenic_regions));
  positions =
      unique_ptr<vector<vector<Gene>>>(new vector<vector<Gene>>(*pi.positions));
  op_max = pi.op_max;
  intergenic_regions_target =
      unique_ptr<vector<IR>>(new vector<IR>(*pi.intergenic_regions_target));
}

Permutation::Permutation(const Genome &g, const Genome &h, bool duplicate) {
  fill(g, h, duplicate);
}

void Permutation::fill(const Genome &g, const Genome &h, bool duplicate) {
  if (duplicate) {
    genes.reset(new vector<Gene>(2 * g.size() - 2));
    intergenic_regions.reset(new vector<IR>(2 * g.size() - 3));
    intergenic_regions_target.reset(new vector<IR>(2 * g.size() - 3));
    op_max = 2 * h.size() - 2;
  } else {
    genes.reset(new vector<Gene>(g.size()));
    intergenic_regions.reset(new vector<IR>(g.size() - 1));
    intergenic_regions_target.reset(new vector<IR>(g.size() - 1));
    op_max = h.size();
  }

  /* Add intergenic regions */
  (*intergenic_regions)[0] = g.get_ir(1);
  for (size_t i = 2; i <= g.size() - 1; i++) {
    if (duplicate) {
      (*intergenic_regions)[2 * i - 3] = 0;
      (*intergenic_regions)[2 * i - 2] = g.get_ir(i);
    } else {
      (*intergenic_regions)[i - 1] = g.get_ir(i);
    }
  }

  /* Add intergenic regions of target genome */
  (*intergenic_regions_target)[0] = h.get_ir(1);
  for (size_t i = 1; i <= h.size() - 1; i++) {
    if (duplicate) {
      (*intergenic_regions_target)[2 * i - 3] = 0;
      (*intergenic_regions_target)[2 * i - 2] = h.get_ir(i);
    } else {
      (*intergenic_regions_target)[i - 1] = h.get_ir(i);
    }
  }

  /* Build vector mapping old to new labels */
  vector<vector<Gene>> labels(h[h.size()] + 1);
  for (size_t i = 0; i < h.size(); i++) {
    labels[abs(h[i + 1])].push_back((h[i + 1] >= 0) ? i : -i);
  }
  for (auto &l : labels) {
    random_shuffle(l.begin(), l.end());
  }

  /* Add genes mapping replicas */
  assert(g[1] >= 0 && g[g.size()] >= 0 && labels[g[1]].size() == 1 &&
         labels[g[g.size()]].size() == 1);
  (*genes)[0] = labels[abs(g[1])].back();
  for (size_t i = 2; i <= g.size() - 1; i++) {
    if (duplicate) {
      if ((g[i] >= 0) == (labels[abs(g[i])].back() >= 0)) {
        (*genes)[2 * i - 3] = 2 * abs(labels[abs(g[i])].back()) - 1;
        (*genes)[2 * i - 2] = 2 * abs(labels[abs(g[i])].back());
      } else {
        (*genes)[2 * i - 3] = 2 * abs(labels[abs(g[i])].back());
        (*genes)[2 * i - 2] = 2 * abs(labels[abs(g[i])].back()) - 1;
      }
    } else {
      (*genes)[i - 1] = labels[abs(g[i])].back();
    }
    labels[abs(g[i])].pop_back();
  }
  if (duplicate) {
    (*genes)[2 * g.size() - 3] = 2*labels[g[g.size()]].back() - 1;
  } else {
    (*genes)[g.size() - 1] = labels[g[g.size()]].back();
  }

  record_positions();
  assert(intergenic_regions->size() == intergenic_regions_target->size());
  assert(intergenic_regions->size() == genes->size() - 1);
}

bool Permutation::is_iota() const {
  bool ok = true;

  for (size_t i = 0; i < size() - 1 && ok; i++) {
    if ((*genes)[i + 1] - (*genes)[i] != 1)
      ok = false;
    if ((*intergenic_regions)[i] != (*intergenic_regions_target)[i])
      ok = false;
  }

  return ok;
}

bool Permutation::breakpoint(Gene ai) const {
  int i = pos(ai)[0];
  return soft_breakpoint(ai) ||
         (*intergenic_regions)[i - 1] !=
             (*intergenic_regions_target)[(*genes)[i - 1]];
}
bool Permutation::hard_breakpoint(Gene ai) const {
  int i = pos(ai)[0];
  return !soft_breakpoint(ai) &&
         (*intergenic_regions)[i - 1] !=
             (*intergenic_regions_target)[(*genes)[i - 1]];
}
bool Permutation::soft_breakpoint(Gene ai) const {
  int i = pos(ai)[0];
  return (*genes)[i] - (*genes)[i - 1] != 1;
}
bool Permutation::overcharged_breakpoint(Gene ai) const {
  int i = pos(ai)[0];
  return !soft_breakpoint(ai) &&
         (*intergenic_regions)[i - 1] >
             (*intergenic_regions_target)[(*genes)[i - 1]];
}

bool Permutation::undercharged_breakpoint(Gene ai) const {
  int i = pos(ai)[0];
  return !soft_breakpoint(ai) &&
         (*intergenic_regions)[i - 1] <
             (*intergenic_regions_target)[(*genes)[i - 1]];
}

void Permutation::serialize(ostream &os) const {
  Genome::serialize(os);
  os << " : [";
  for (size_t i = 1; i < size() - 1; ++i) {
    os << get_ir_target(i) << ", ";
  }
  os << get_ir_target(size() - 1) << "]";
  /* for (size_t i = 1; i <= size(); ++i) { */
  /*   os << (*this)[i] << " "; */
  /* } */
  /* os << endl; */
  /* for (size_t i = 1; i < size(); ++i) { */
  /*   os << get_ir(i) << " "; */
  /* } */
  /* os << endl; */
  /* for (size_t i = 1; i <= size(); ++i) { */
  /*   os << i-1 << " "; */
  /* } */
  /* os << endl; */
  /* for (size_t i = 1; i < size(); ++i) { */
  /*   os << get_ir_target(i) << " "; */
  /* } */
}
