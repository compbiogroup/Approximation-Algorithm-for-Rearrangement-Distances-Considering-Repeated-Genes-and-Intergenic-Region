#include "transposition3.hpp"
#include "../misc/permutation.hpp"
#include "aux.hpp"
#include <algorithm>
#include <cassert>
#include <set>

int Transposition3::estimate_distance(Permutation pi) {
  int i, j, k, dist = 0;
  set<IR> breaks;
  bool single_excess =
      false; // When set to true we know that there is a single breakpoint with
             // excess of nucleotides or there is no soft breakpoints

  for (size_t b = 1; b <= pi.size() - 1; b++) {
    if (pi.breakpoint(pi[b], false))
      breaks.insert(pi[b]);
  }
  size_t breaks_size = breaks.size() + 1;

  while (!breaks.empty()) {

    assert(breaks.size() < breaks_size);
    /* if (breaks.size() >= breaks_size) { */
    /* return -1; */
    /* } */
    breaks_size = breaks.size();

    // If there is two or more overcharged breakpoints
    bool ok = false;
    vector<Gene> over;
    Gene other = -1;

    for (auto b : breaks) {
      if (over.size() < 2 && pi.overcharged_breakpoint(b, false)) {
        over.push_back(b);
      } else {
        other = b;
      }
      if (over.size() >= 2 && other != -1) {
        ok = true;
        break;
      }
    }

    if (ok) {
      vector<pair<int, IR>> idxs_ir;
      IR old_sum =
          pi.get_ir(pi.pos(over[0])[0]) + pi.get_ir(pi.pos(over[1])[0]);
      IR new_sum =
          pi.get_ir_target(over[0] + 1) + pi.get_ir_target(over[1] + 1);
      idxs_ir.push_back(
          pair<int, IR>(pi.pos(over[0])[0], pi.get_ir_target(over[0] + 1)));
      idxs_ir.push_back(
          pair<int, IR>(pi.pos(over[1])[0], pi.get_ir_target(over[1] + 1)));
      idxs_ir.push_back(pair<int, IR>(
          pi.pos(other)[0], pi.get_ir(pi.pos(other)[0]) + old_sum - new_sum));
      sort(idxs_ir.begin(), idxs_ir.end());
      transposition_move(pi, idxs_ir[0].first, idxs_ir[1].first,
                         idxs_ir[2].first, idxs_ir[0].second, idxs_ir[1].second,
                         idxs_ir[2].second);
      dist += 2;
      for (auto b : over) {
        if (not pi.breakpoint(b, false))
          breaks.erase(b);
      }
      if (not pi.breakpoint(other, false))
        breaks.erase(other);
      continue;
    }

    // If there exists a soft breakpoint with excess of nucleotides
    int a, c;
    for (auto b : breaks) {
      if (pi.soft_breakpoint(b, false) &&
          pi.get_ir(pi.pos(b)[0]) >= pi.get_ir_target(b + 1)) {
        i = pi.pos(b)[0];
        a = b;
        ok = true;
        break;
      }
    }
    if (ok) {
      j = pi.pos(pi[i] + 1)[0];
      c = pi[j - 1];
      if (j > i) {
        for (auto b : breaks) {
          k = pi.pos(b)[0];
          if (k < i) {
            // If single_excess is true, there exists at least 4 soft
            // breakpoints and the region between S_{k} and S_{i+1} would become
            // an overcharged breakpoints then we change the transposition
            if (single_excess && breaks.size() >= 4 && pi[i + 1] - pi[k] == 1) {
              int d, l = -1;
              for (auto d_ : breaks) {
                if (d_ != a && d_ != b && d_ != c) {
                  d = d_;
                  l = pi.pos(d)[0];
                  break;
                }
              }
              assert(l != -1);
              if (l < k) {
                pi.transposition(l + 1, i + 1, j, 0,
                                 pi.get_ir_target(pi[i] + 1), pi.get_ir(j - 1));
              } else if (l >= j) {
                pi.transposition(i + 1, j, l + 1, pi.get_ir_target(pi[i] + 1),
                                 pi.get_ir(j - 1), 0);
              } else if (l < i) { // and l > k
                pi.transposition(l + 1, i + 1, j, 0,
                                 pi.get_ir_target(pi[i] + 1), pi.get_ir(j - 1));
              } else { // l > i and l < j-1
                pi.transposition(k + 1, i + 1, l + 1, 0,
                                 pi.get_ir(i) - pi.get_ir_target(pi[k] + 1), 0);
              }
              if (not pi.breakpoint(d, false))
                breaks.erase(d);
            } else {
              pi.transposition(k + 1, i + 1, j, 0, pi.get_ir_target(pi[i] + 1),
                               pi.get_ir(j - 1));
            }
            if (not pi.breakpoint(a, false))
              breaks.erase(a);
            if (not pi.breakpoint(b, false))
              breaks.erase(b);
            if (not pi.breakpoint(c, false))
              breaks.erase(c);
            break;
          } else if (k >= j) {
            // If single_excess is true, there exists at least 4 soft
            // breakpoints and the region between S_{k} and S_{i+1} would become
            // an overcharged breakpoints then we change the transposition
            if (single_excess && breaks.size() >= 4 && pi[i + 1] - pi[k] == 1) {
              int d, l = -1;
              for (auto d_ : breaks) {
                if (d_ != a && d_ != b && d_ != c) {
                  d = d_;
                  l = pi.pos(d)[0];
                  break;
                }
              }
              assert(l != -1);
              if (l < i) {
                pi.transposition(l + 1, i + 1, j, 0,
                                 pi.get_ir_target(pi[i] + 1), pi.get_ir(j - 1));
              } else if (l >= j) {
                pi.transposition(i + 1, j, l + 1, pi.get_ir_target(pi[i] + 1),
                                 pi.get_ir(j - 1), 0);
              } else { // l < j - 1 and l > i
                pi.transposition(i + 1, l + 1, k + 1,
                                 pi.get_ir(i) - pi.get_ir_target(pi[k] + 1), 0,
                                 0);
              }
              if (not pi.breakpoint(d, false))
                breaks.erase(d);
            } else {
              pi.transposition(i + 1, j, k + 1, pi.get_ir_target(pi[i] + 1),
                               pi.get_ir(j - 1), 0);
            }
            if (not pi.breakpoint(a, false))
              breaks.erase(a);
            if (not pi.breakpoint(b, false))
              breaks.erase(b);
            if (not pi.breakpoint(c, false))
              breaks.erase(c);
            break;
          }
        }
      } else {
        for (auto b : breaks) {
          k = pi.pos(b)[0];
          if (j <= k && k < i) {
            // If single_excess is true, there exists at least 4 soft
            // breakpoints and the region between S_{k} and S_{i+1} would become
            // an overcharged breakpoints then we change the transposition
            if (single_excess && breaks.size() >= 4 && pi[i + 1] - pi[k] == 1) {
              int d, l = -1;
              for (auto d_ : breaks) {
                if (d_ != a && d_ != b && d_ != c) {
                  d = d_;
                  l = pi.pos(d)[0];
                  break;
                }
              }
              assert(l != -1);
              if (l < j - 1) {
                pi.transposition(l + 1, k + 1, i + 1, 0, 0,
                                 pi.get_ir(i) - pi.get_ir_target(pi[k] + 1));
              } else if (l > i) {
                pi.transposition(j, i + 1, l + 1, pi.get_ir(j - 1),
                                 pi.get_ir_target(pi[i] + 1), 0);
              } else { // l >= j and l < i
                pi.transposition(j, l + 1, i + 1, pi.get_ir(j - 1), 0,
                                 pi.get_ir_target(pi[i] + 1));
              }
              if (not pi.breakpoint(d, false))
                breaks.erase(d);
            } else {
              pi.transposition(j, k + 1, i + 1, pi.get_ir(j - 1), 0,
                               pi.get_ir_target(pi[i] + 1));
            }
            if (not pi.breakpoint(a, false))
              breaks.erase(a);
            if (not pi.breakpoint(b, false))
              breaks.erase(b);
            if (not pi.breakpoint(c, false))
              breaks.erase(c);
            break;
          }
        }
      }
      dist += 1;
      continue;
    }

    // There is only one overcharged breakpoint and no soft breakpoint has extra
    // nucleotides
    assert(!over.empty());
    Gene under = -1;
    for (auto b : breaks) {
      if (pi.undercharged_breakpoint(b, false)) {
        under = b;
        break;
      }
    }
    if (under != -1) {
      other = -1;
      for (auto b : breaks) {
        if (b != under && b != over[0]) {
          other = b;
          break;
        }
      }
      IR old_sum = pi.get_ir(pi.pos(under)[0]) + pi.get_ir(pi.pos(over[0])[0]);
      IR new_sum = pi.get_ir_target(under + 1) + pi.get_ir_target(over[0] + 1);
      if (other == -1) {
        assert(old_sum == new_sum);
        for (size_t b = 1; b <= pi.size() - 1; b++) {
          if (!pi.breakpoint(pi[b], false))
            other = pi[b];
        }
      } else {
        assert(old_sum >= new_sum);
      }
      vector<pair<int, IR>> idxs_ir;
      idxs_ir.push_back(
          pair<int, IR>(pi.pos(over[0])[0], pi.get_ir_target(over[0] + 1)));
      idxs_ir.push_back(
          pair<int, IR>(pi.pos(under)[0], pi.get_ir_target(under + 1)));
      idxs_ir.push_back(pair<int, IR>(
          pi.pos(other)[0], pi.get_ir(pi.pos(other)[0]) + old_sum - new_sum));
      sort(idxs_ir.begin(), idxs_ir.end());
      transposition_move(pi, idxs_ir[0].first, idxs_ir[1].first,
                         idxs_ir[2].first, idxs_ir[0].second, idxs_ir[1].second,
                         idxs_ir[2].second);
      dist += 2;
      if (not pi.breakpoint(over[0], false))
        breaks.erase(over[0]);
      if (not pi.breakpoint(under, false))
        breaks.erase(under);
      if (not pi.breakpoint(other, false))
        breaks.erase(other);
    } else {
      vector<Gene> softs;
      for (auto b : breaks) {
        if (b != over[0]) {
          softs.push_back(b);
        }
        if (softs.size() >= 2)
          break;
      }
      vector<pair<int, IR>> idxs_ir;
      IR old_sum =
          pi.get_ir(pi.pos(softs[0])[0]) + pi.get_ir(pi.pos(over[0])[0]);
      IR new_sum =
          pi.get_ir_target(softs[0] + 1) + pi.get_ir_target(over[0] + 1);
      idxs_ir.push_back(
          pair<int, IR>(pi.pos(over[0])[0], pi.get_ir_target(over[0] + 1)));
      idxs_ir.push_back(
          pair<int, IR>(pi.pos(softs[0])[0], pi.get_ir_target(softs[0] + 1)));
      idxs_ir.push_back(
          pair<int, IR>(pi.pos(softs[1])[0],
                        pi.get_ir(pi.pos(softs[1])[0]) + old_sum - new_sum));
      sort(idxs_ir.begin(), idxs_ir.end());
      transposition_move(pi, idxs_ir[0].first, idxs_ir[1].first,
                         idxs_ir[2].first, idxs_ir[0].second, idxs_ir[1].second,
                         idxs_ir[2].second);
      dist += 2;
      if (not pi.breakpoint(over[0], false))
        breaks.erase(over[0]);
      if (not pi.breakpoint(softs[0], false))
        breaks.erase(softs[0]);
      if (not pi.breakpoint(softs[1], false))
        breaks.erase(softs[1]);
      single_excess = true;
    }
  }

  assert(pi.is_iota());
  return dist;
}
