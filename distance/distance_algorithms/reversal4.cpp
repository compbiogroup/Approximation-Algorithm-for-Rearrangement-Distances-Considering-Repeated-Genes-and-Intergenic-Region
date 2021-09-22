#include "reversal4.hpp"
#include "../misc/permutation.hpp"
#include "aux.hpp"
#include <cassert>
#include <limits>
#include <vector>

int Reversal4::estimate_distance(Permutation pi) {
  int i, j, dist = 0, step, case_x = 0, a_pos, b_pos;
  vector<IR> breaks;
  size_t breaks_size = numeric_limits<int>::max();

  while (true) {

    breaks.clear();
    for (size_t b = 1; b <= pi.size() - 1; b++) {
      if (pi.breakpoint(pi[b], true)) {
        breaks.push_back(pi[b]);
      }
    }

    if (breaks.empty())
      break;

    assert(breaks.size() < breaks_size);
    breaks_size = breaks.size();

    case_x = 0;
    for (auto a : breaks) {
      for (auto b : breaks) {
        a_pos = pi.pos(a)[0];
        b_pos = pi.pos(b)[0];
        if (a_pos < b_pos) {
          step = is_connected(pi, a_pos, b_pos);
          if (step != 0 && (case_x == 0 || case_x > step)) {
            case_x = step;
            i = a_pos;
            j = b_pos;
          }
        }
      }
    }

    if (case_x == 1) { // (i, j) / (i + 1, j + 1)

      dist += case_i(pi, i, j);

    } else if (case_x == 2) { // (i, j + 1)

      for (auto c : breaks) {
        int k = pi.pos(c)[0];
        if (k < i) {
          pi.reversal(k + 1, i, 0, pi.get_ir(i));
          dist += 1 + case_i(pi, k, j);
          break;
        }
        if (k > j) {
          pi.reversal(j + 1, k, 0, pi.get_ir(k));
          dist += 1 + case_i(pi, i, k);
          break;
        }
      }

    } else if (case_x == 3) { // (i + 1, j)

      for (auto c : breaks) {
        int k = pi.pos(c)[0];
        if (k > i && k < j) {
          pi.reversal(i + 1, k, 0, pi.get_ir(k));
          dist += 1 + case_i(pi, k, j);
          break;
        }
      }

    } else if (case_x == 4) { // (i, i + 1) / (j, j + 1)

      pi.reversal(i + 1, j, 0, 0);
      dist += 1 + case_i(pi, i, j);
    }
  }

  assert(pi.is_iota());
  return dist;
}
