#include <algorithm>
#include <vector>

#include "Transposition4/transposition4.hpp"
#include "external/external.hpp"
#include "misc/genome.hpp"
#include "misc/permutation.hpp"
#include "quickcheck/quickcheck/quickcheck.hh"
using namespace quickcheck;

void generate(size_t n, Permutation &perm) {
  string sg, ig, sh, ih;
  int ir;
  int sum = 0;
  bool add_sign, positive;

  generate(n, add_sign);
  sg.append("0");
  sg.append(" ");
  sh.append("0");
  sh.append(" ");
  ir = generateInRange(50, 100);
  ig.append(to_string(ir));
  ig.append(" ");
  sum += ir;
  for (size_t i = 0; i < n; i++) {
    ir = generateInRange(0, 100);
    generate(n, positive);
    if (!add_sign) {
      sg.append("1");
    } else if (positive) {
      sg.append("+1");
    } else {
      sg.append("-1");
    }
    sg.append(" ");
    ig.append(to_string(ir));
    ig.append(" ");
    generate(n, positive);
    if (!add_sign) {
      sh.append("1");
    } else if (positive) {
      sh.append("+1");
    } else {
      sh.append("-1");
    }
    sh.append(" ");
    sum += ir;
  }
  sg.append("2");
  sg.append(" ");
  sh.append("2");
  sh.append(" ");

  vector<int> ir_list;
  ir_list.push_back(0);
  for (size_t i = 0; i < n; i++) {
    ir_list.push_back(generateInRange(0, sum));
  }
  ir_list.push_back(sum);
  sort(ir_list.begin(), ir_list.end());
  for (size_t i = 1; i <= n + 1; i++) {
    ir = ir_list[i] - ir_list[i - 1];
    ih.append(to_string(ir));
    ih.append(" ");
  }

  Genome g(sg, ig, false);
  Genome h(sh, ih, false);
  perm.fill(g, h, add_sign);
}

class PBounds4T : public Property<Permutation> {
  bool holdsFor(const Permutation &perm) {
    Transposition4 trans_alg;

    int breaks_count = 0;
    for (size_t b = 1; b <= perm.size() - 1; b++) {
      if (perm.breakpoint(perm[b]))
        breaks_count++;
    }

    int dist = trans_alg.estimate_distance(perm);

    bool sucess = true;
    sucess = sucess && dist >= (breaks_count / 3);
    sucess = sucess && dist <= (4 * breaks_count / 3);
    if (!sucess) {
      cout << "pi: " << perm << endl;
      cout << "dist: " << dist << endl;
      cout << "b(pi):" << breaks_count << endl;
    } else {
      cout << ".";
      cout.flush();
    }
    return sucess;
  }
};


int main() {
  // set seed
  srand(time(nullptr));

  check<PBounds4T>(
      "upper and lower bounds hold for factor 4 algorithm for transposition");
  return 0;
}
