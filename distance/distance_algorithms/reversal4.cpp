#include "reversal4.hpp"
#include "../misc/permutation.hpp"
#include <cassert>
#include <limits>
#include <set>

int is_connected(const Permutation &pi, int i, int j) {
  // return values:
  // 0 - not connected
  // 1 - connected (case i)
  // 2 - connected (case ii)
  // 3 - connected (case iii)
  // 4 - connected (case iv)
  int nucleotides = pi.get_ir(i) + pi.get_ir(j);

  if (abs(pi[i] - pi[j]) == 1 && !pi.is_block(pi[i], pi[j], true) &&
      pi.get_ir_target(min(pi[i], pi[j]) + 1) <= nucleotides) {
    return 1;
  }
  if (abs(pi[i + 1] - pi[j + 1]) == 1 &&
      !pi.is_block(pi[i + 1], pi[j + 1], true) &&
      pi.get_ir_target(min(pi[i + 1], pi[j + 1]) + 1) <= nucleotides) {
    return 1;
  }
  if (abs(pi[i] - pi[j + 1]) == 1 && !pi.is_block(pi[i], pi[j + 1], true) &&
      pi.get_ir_target(min(pi[i], pi[j + 1]) + 1) <= nucleotides) {
    return 2;
  }
  if (abs(pi[i + 1] - pi[j]) == 1 && !pi.is_block(pi[i + 1], pi[j], true) &&
      pi.get_ir_target(min(pi[i + 1], pi[j]) + 1) <= nucleotides) {
    return 3;
  }
  if (abs(pi[i] - pi[i + 1]) == 1 && !pi.is_block(pi[i], pi[i + 1], true) &&
      pi.get_ir_target(min(pi[i], pi[i + 1]) + 1) <= nucleotides) {
    return 4;
  }
  if (abs(pi[j] - pi[j + 1]) == 1 && !pi.is_block(pi[j], pi[j + 1], true) &&
      pi.get_ir_target(min(pi[j], pi[j + 1]) + 1) <= nucleotides) {
    return 4;
  }
  return 0;
}

int case_1(Permutation &pi, int i, int j) {
  int goal, x, y;
  if (abs(pi[i] - pi[j]) == 1) {
    goal = pi.get_ir_target(min(pi[i], pi[j]) + 1);
    if (pi.get_ir(i) + pi.get_ir(j) >= goal) {
      if (pi.get_ir(i) >= goal) {
        x = goal;
        y = 0;
        pi.reversal(i + 1, j, x, y);
      } else {
        x = pi.get_ir(i);
        y = goal - x;
        pi.reversal(i + 1, j, x, y);
      }
      return 1;
    }
  }
  if (abs(pi[i + 1] - pi[j + 1]) == 1) {
    goal = pi.get_ir_target(min(pi[i + 1], pi[j + 1]) + 1);
    if (pi.get_ir(i) + pi.get_ir(j) >= goal) {
      if (pi.get_ir(j) >= goal) {
        x = pi.get_ir(i);
        y = pi.get_ir(j) - goal;
        pi.reversal(i + 1, j, x, y);
      } else {
        x = pi.get_ir(i) - (goal - pi.get_ir(j));
        y = 0;
        pi.reversal(i + 1, j, x, y);
      }
      return 1;
    }
  }
  return 0;
}

int Reversal4::estimate_distance(Permutation pi) {
  int i, j, dist = 0, step, case_x = 0, a_pos, b_pos;
  set<IR> breaks;
  size_t breaks_size = numeric_limits<int>::max();

  while (true) {

    breaks.clear();
    for (size_t b = 1; b <= pi.size() - 1; b++) {
      if (pi.breakpoint(pi[b], true)) {
        breaks.insert(pi[b]);
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

      dist += case_1(pi, i, j);

    } else if (case_x == 2) { // (i, j + 1)

      for (auto c : breaks) {
        int k = pi.pos(c)[0];
        if (k < i) {
          pi.reversal(k + 1, i, 0, pi.get_ir(i));
          dist += 1 + case_1(pi, k, j);
          break;
        }
        if (k > j) {
          pi.reversal(j + 1, k, 0, pi.get_ir(k));
          dist += 1 + case_1(pi, i, k);
          break;
        }
      }

    } else if (case_x == 3) { // (i + 1, j)

      for (auto c : breaks) {
        int k = pi.pos(c)[0];
        if (k > i && k < j) {
          pi.reversal(i + 1, k, 0, pi.get_ir(k));
          dist += 1 + case_1(pi, k, j);
          break;
        }
      }

    } else if (case_x == 4) { // (i, i + 1) / (j, j + 1)

      pi.reversal(i + 1, j, 0, 0);
      dist += 1 + case_1(pi, i, j);
    }
  }

  assert(pi.is_iota());
  return dist;
}
