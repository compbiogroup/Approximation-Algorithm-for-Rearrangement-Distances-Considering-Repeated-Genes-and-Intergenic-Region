#include "io.hpp"
#include "genome.hpp"
#include "../cycle/cycles.hpp"
#include <iostream>
#include <sstream>

vector<string> *read_lines(istream &is) {
  auto input_lines = new vector<string>();

  for (string line; getline(is, line);) {
    if (line.at(0) != '#') {
      input_lines->push_back(line);
    }
  }

  return input_lines;
}

InputData input(string &line1, string &line2, string &line3, string &line4, bool extend) {
  string str;
  InputData data;

  data.g = new Genome(line1, line2, extend);
  data.h = new Genome(line3, line4, extend);
  if (data.g->size() != data.h->size()) {
    throw invalid_argument("Size of genomes differ");
  }

  return data;
}

void output(ostream &os, int dist, double time) {
  os << "Dist: " << dist;
  os.precision(5);
  os << fixed;
  os << ", Wall Time: " << time << "s" << endl;
}

void output(ostream &os, const CycleGraph &cyc_dec, double time) {
  os << "Dec: " << cyc_dec;
  os.precision(5);
  os << fixed;
  os << ", Wall Time: " << time << "s" << endl;
}

void output(ostream &os, PermsIrs perms_irs) {
    for (size_t i = 0; i < perms_irs.n - 1; i++) {
        os << perms_irs.s[i] << " ";
    }
    os << perms_irs.s[perms_irs.n - 1] << endl;
    for (size_t i = 0; i < perms_irs.n - 2; i++) {
        os << perms_irs.s_ir[i] << " ";
    }
    os << perms_irs.s_ir[perms_irs.n - 2] << endl;
    for (size_t i = 0; i < perms_irs.n - 1; i++) {
        os << perms_irs.p[i] << " ";
    }
    os << perms_irs.p[perms_irs.n - 1] << endl;
    for (size_t i = 0; i < perms_irs.n - 2; i++) {
        os << perms_irs.p_ir[i] << " ";
    }
    os << perms_irs.p_ir[perms_irs.n - 2] << endl;
}
