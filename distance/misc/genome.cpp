#include "genome.hpp"
#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <sstream>

Genome::Genome(string str_g, string str_i, bool extend) {
  string token;
  op_max = -1;
  genes.reset(new vector<Gene>());
  intergenic_regions.reset(new vector<IR>());

  /* Read each element of the string. */
  if (extend)
    genes->push_back(0);
  stringstream ssg(str_g);
  while (getline(ssg, token, ' ')) {
    Gene a = stoi(token);
    genes->push_back(a);
    if (abs(genes->back()) > op_max)
      op_max = abs(genes->back());
  }
  if (extend) {
    genes->push_back(op_max + 1);
    op_max++;
  }

  /* Read each intergenic regions. */
  stringstream ssi(str_i);
  while (getline(ssi, token, ' ')) {
    IR r = stoi(token);
    intergenic_regions->push_back(r);
  }

  record_positions();

  assert(intergenic_regions->size() == genes->size() - 1);
}

void Genome::record_positions() {
  /* Record list of positions for each label. */
  positions.reset(new vector<vector<Gene>>(op_max + 1));
  for (size_t i = 0; i < genes->size(); ++i) {
    (*positions)[abs((*genes)[i])].push_back(i + 1);
  }
}

int Genome::occ_max() const {
  vector<int> count(this->positions->size(), 0);
  for (size_t i = 0; i < this->size(); ++i) {
    count[abs((*this)[i])]++;
  }
  return *max_element(count.begin(), count.end());
}

void Genome::deletion(Gene i, Gene j, IR x) {
  assert(2 <= i);
  assert(i < j);
  assert(j <= int(size()));
  assert(0 <= x &&
         x <= (*intergenic_regions)[i - 2] + (*intergenic_regions)[j - 2]);

  int k = i - 1, l = j - 1;
  (*intergenic_regions)[i - 2] = x;
  for (; l < int(size()) - 1; k++, l++) {
    (*genes)[k] = (*genes)[l];
    (*intergenic_regions)[k] = (*intergenic_regions)[l];
  }
  (*genes)[k] = (*genes)[l];
  genes->resize(k+1);
  intergenic_regions->resize(k);
  record_positions();
}

void Genome::reversal(Gene i, Gene j, IR x, IR y) {

  assert(2 <= i);
  assert(i <= j);
  assert(j <= int(size()));
  assert(0 <= x && x <= (*intergenic_regions)[i - 2]);
  assert(0 <= y && y <= (*intergenic_regions)[j - 1]);

  int x_rest = (*intergenic_regions)[i - 2] - x;
  int y_rest = (*intergenic_regions)[j - 1] - y;

  for (int k = 0; k < (j - i + 1) / 2; k++) {
    swap((*genes)[i + k - 1], (*genes)[j - k - 1]);
    swap((*intergenic_regions)[i + k - 1], (*intergenic_regions)[j - k - 2]);
  }

  (*intergenic_regions)[i - 2] = x + y;
  (*intergenic_regions)[j - 1] = x_rest + y_rest;
  record_positions();
}

void Genome::transposition(Gene i, Gene j, Gene k, IR x, IR y, IR z) {

  assert(2 <= i);
  assert(i < j);
  assert(j < k);
  assert(k <= int(size()));
  assert(0 <= x && x <= (*intergenic_regions)[i - 2]);
  assert(0 <= y && y <= (*intergenic_regions)[j - 2]);
  assert(0 <= z && z <= (*intergenic_regions)[k - 2]);

  int x_rest = (*intergenic_regions)[i - 2] - x;
  int y_rest = (*intergenic_regions)[j - 2] - y;
  int z_rest = (*intergenic_regions)[k - 2] - z;

  vector<Gene> aux1;
  vector<IR> aux2;
  for (int t = i - 1; t <= j - 2; t++) {
    aux1.push_back((*genes)[t]);
    aux2.push_back((*intergenic_regions)[t]);
  }
  int t = i - 1;
  for (int r = j - 1; r <= k - 2; r++) {
    (*genes)[t] = (*genes)[r];
    (*intergenic_regions)[t] = (*intergenic_regions)[r];
    t++;
  }
  for (size_t r = 0; r < aux1.size(); r++) {
    (*genes)[t] = aux1[r];
    (*intergenic_regions)[t] = aux2[r];
    t++;
  }

  (*intergenic_regions)[i - 2] = x + y_rest;
  (*intergenic_regions)[i + k - j - 2] = z + x_rest;
  (*intergenic_regions)[k - 2] = y + z_rest;
  record_positions();
}

void Genome::serialize(ostream &os) const {
  for (size_t i = 1; i <= size() - 1; ++i) {
    os << "(" << (*this)[i] << ") - " << get_ir(i) << " - ";
  }
  os << "(" << (*this)[size()] << ")";
}

ostream &operator<<(ostream &os, const Genome &g) {
  g.serialize(os);
  return os;
}
