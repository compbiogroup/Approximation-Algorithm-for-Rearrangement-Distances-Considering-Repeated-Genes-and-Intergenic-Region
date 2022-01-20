#include "reversal_transposition45.hpp"
#include "../misc/permutation.hpp"
#include "aux.hpp"
#include "../misc/io.hpp"
#include <cassert>
#include <limits>
#include <vector>
// (i, j + 1)
int case_ii(Permutation &pi, vector<int> &breaks,
            vector<int> &overcharged_breaks, int i, int j) {
  int k, x, y, z, goal = pi.get_ir_target(min(pi[i], pi[j + 1]) + 1);
  for (auto k_it = breaks.begin(); k_it != overcharged_breaks.end(); k_it++) {
    if (k_it == breaks.end())
      k_it = overcharged_breaks.begin();
    k = pi.pos(*k_it)[0];
    if (k < i) {
      if (pi.get_ir(j) >= goal) {
        x = 0;
        y = 0;
        z = pi.get_ir(j) - goal;
        pi.transposition(k + 1, i + 1, j + 1, x, y, z);
        return 1;
      } else {
        x = 0;
        y = goal - pi.get_ir(j);
        z = 0;
        pi.transposition(k + 1, i + 1, j + 1, x, y, z);
        return 1;
      }
    } else if (k > j) {
      if (pi.get_ir(i) >= goal) {
        x = goal;
        y = pi.get_ir(j);
        z = 0;
        pi.transposition(i + 1, j + 1, k + 1, x, y, z);
        return 1;
      } else {
        x = pi.get_ir(i);
        y = pi.get_ir(j) - (goal - pi.get_ir(i));
        z = 0;
        pi.transposition(i + 1, j + 1, k + 1, x, y, z);
        return 1;
      }
    }
  }
  return 0;
}

// (i + 1, j)
int case_iii(Permutation &pi, vector<int> &breaks,
             vector<int> &overcharged_breaks, int i, int j) {
  int k, x, y, z, goal = pi.get_ir_target(min(pi[i + 1], pi[j]) + 1);
  for (auto k_it = breaks.begin(); k_it != overcharged_breaks.end(); k_it++) {
    if (k_it == breaks.end())
      k_it = overcharged_breaks.begin();
    k = pi.pos(*k_it)[0];
    if (k > i and k < j) {
      if (pi.get_ir(j) >= goal) {
        x = pi.get_ir(i);
        y = 0;
        z = goal;
        pi.transposition(i + 1, k + 1, j + 1, x, y, z);
        return 1;
      } else {
        x = pi.get_ir(i) - (goal - pi.get_ir(j));
        y = 0;
        z = pi.get_ir(j);
        pi.transposition(i + 1, k + 1, j + 1, x, y, z);
        return 1;
      }
    }
  }
  return 0;
}

// (i, i + 1) / (j, j + 1)
int case_iv(Permutation &pi, int i, int j) {
  pi.reversal(i + 1, j, 0, 0);
  return 1 + case_i(pi, i, j);
}

int two_overcharged(Permutation &pi, vector<int> breaks, vector<int> overcharged_breaks) {
    int i,j,k,x,y,z,a,b,c,b0,ob0,ob1;

    b0 = pi.pos(breaks[0])[0];
    ob0 = pi.pos(overcharged_breaks[0])[0];
    ob1 = pi.pos(overcharged_breaks[1])[0];
    if (b0 < ob0) {

      i = b0;
      j = ob0;
      k = ob1;

      a = pi.get_ir(i);
      b = pi.get_ir(j);
      c = pi.get_ir(k);

      y = pi.get_ir_target(min(pi[ob0], pi[ob0 + 1]) + 1);
      z = pi.get_ir_target(min(pi[ob1], pi[ob1 + 1]) + 1);
      x = a + (b - y) + (c - z);

      transposition_move(pi, i, j, k, x, y, z);
      return 2;

    } else if (b0 > ob1) {

      i = ob0;
      j = ob1;
      k = b0;

      a = pi.get_ir(i);
      b = pi.get_ir(j);
      c = pi.get_ir(k);

      x = pi.get_ir_target(min(pi[ob0], pi[ob0 + 1]) + 1);
      y = pi.get_ir_target(min(pi[ob1], pi[ob1 + 1]) + 1);
      z = c + (a - x) + (b - y);

      transposition_move(pi, i, j, k, x, y, z);
      return 2;

    } else {

      i = ob0;
      j = b0;
      k = ob1;

      a = pi.get_ir(i);
      b = pi.get_ir(j);
      c = pi.get_ir(k);

      x = pi.get_ir_target(min(pi[ob0], pi[ob0 + 1]) + 1);
      z = pi.get_ir_target(min(pi[ob1], pi[ob1 + 1]) + 1);
      y = b + (a - x) + (c - z);

      transposition_move(pi, i, j, k, x, y, z);
      return 2;
    }
}

int one_overcharged(Permutation &pi, vector<int> breaks, vector<int> overcharged_breaks) {
    int i,j,k,a,b,c,x,y,z,b0,b1,ob0;

    b0 = pi.pos(breaks[0])[0];
    ob0 = pi.pos(overcharged_breaks[0])[0];
    if (pi.undercharged_breakpoint(breaks[0], true)) {

      a = pi.get_ir(ob0);
      b = pi.get_ir(b0);
      x = pi.get_ir_target(min(pi[ob0], pi[ob0 + 1]) + 1);
      y = pi.get_ir_target(min(pi[b0], pi[b0 + 1]) + 1);

      if ((a + b) > (x + y)) {

        b1 = pi.pos(breaks[1])[0];
        if (ob0 < b0) {

          i = ob0;
          j = b0;
          k = b1;

          a = pi.get_ir(i);
          b = pi.get_ir(j);
          c = pi.get_ir(k);

          x = pi.get_ir_target(min(pi[ob0], pi[ob0 + 1]) + 1);
          y = pi.get_ir_target(min(pi[b0], pi[b0 + 1]) + 1);
          z = (a + b + c) - (x + y);

          transposition_move(pi, i, j, k, x, y, z);
          return 2;

        }  else if (ob0 > b1) {

          i = b0;
          j = b1;
          k = ob0;

          a = pi.get_ir(i);
          b = pi.get_ir(j);
          c = pi.get_ir(k);

          x = pi.get_ir_target(min(pi[b0], pi[b0 + 1]) + 1);
          z = pi.get_ir_target(min(pi[ob0], pi[ob0 + 1]) + 1);
          y = (a + b + c) - (x + z);

          transposition_move(pi, i, j, k, x, y, z);
          return 2;

        } else {

          i = b0;
          j = ob0;
          k = b1;

          a = pi.get_ir(i);
          b = pi.get_ir(j);
          c = pi.get_ir(k);

          x = pi.get_ir_target(min(pi[b0], pi[b0 + 1]) + 1);
          y = pi.get_ir_target(min(pi[ob0], pi[ob0 + 1]) + 1);
          z = (a + b + c) - (x + y);

          transposition_move(pi, i, j, k, x, y, z);
          return 2;
        }
      } else { 

        if (b0 > ob0) {
          return case_iv(pi, ob0, b0);
        } else {
          return case_iv(pi, b0, ob0);
        }
      }
    } else {

      if (b0 > ob0) {
        return case_iv(pi, ob0, b0);
      } else { 
        return case_iv(pi, b0, ob0);
      }
    }
}

int ReversalTransposition45::estimate_distance(Permutation pi) {
  int i, j, dist = 0, step = 0, stage_1 = true, case_x = 0, a_pos, b_pos;
  vector<IR> breaks, overcharged_breaks;
  size_t breaks_size = numeric_limits<int>::max();
  bool swap_stage = false;

  while (true) {

    breaks.clear();
    overcharged_breaks.clear();
    for (size_t b = 1; b <= pi.size() - 1; b++) {
      if (pi.overcharged_breakpoint(pi[b], true)) {
        overcharged_breaks.push_back(pi[b]);
      } else if (pi.breakpoint(pi[b], true)) {
        breaks.push_back(pi[b]);
      }
    }
    /* cout << pi << endl; */
    /* cout << endl << "b: "; */
    /* print_vec(cout, breaks); */
    /* cout << endl << "bo: "; */
    /* print_vec(cout, overcharged_breaks); */
    /* cout << endl; */

    if (overcharged_breaks.empty() && breaks.empty())
      break;

    if (swap_stage) {
        stage_1 = not stage_1;
    } else {
        assert(breaks.size() + overcharged_breaks.size() < breaks_size);
        breaks_size = breaks.size() + overcharged_breaks.size();
    }

    swap_stage = true;
    if (stage_1) {
      if (overcharged_breaks.size() > 1) {
        dist += two_overcharged(pi, breaks, overcharged_breaks);
        swap_stage = false;
      } else if (overcharged_breaks.size() == 1) {
        dist += one_overcharged(pi, breaks, overcharged_breaks);
        swap_stage = false;
      }
    } else {
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

      if (case_x == 1) {
        dist += case_i(pi, i, j);
        swap_stage = false;
      } else if (case_x == 2) {
        dist += case_ii(pi, breaks, overcharged_breaks, i, j);
        swap_stage = false;
      } else if (case_x == 3) {
        dist += case_iii(pi, breaks, overcharged_breaks, i, j);
        swap_stage = false;
      }
    }
  }

  assert(pi.is_iota());
  return dist;
}
