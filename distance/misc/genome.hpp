#pragma once

#include <iostream>
#include <memory>
#include <vector>
using namespace std;

typedef int Gene;
typedef int IR;

/* String representing a unichromossomal and linear genome */
class Genome {
protected:
  unique_ptr<vector<Gene>> genes;
  unique_ptr<vector<IR>> intergenic_regions;
  unique_ptr<vector<vector<Gene>>> positions;
  Gene op_max;
  Genome(){};
  void record_positions();

public:
  Genome(string str_g, string str_i, bool extend);
  size_t size() const { return genes->size(); }
  /* Only read access to the elements. */
  Gene operator[](int i) const { return (*genes)[i - 1]; }
  /* Only read access to intergenic regions. */
  IR get_ir(int i) const { return (*intergenic_regions)[i - 1]; }
  /* Get maximum occurrence */
  int occ_max() const;
  /* Get positions of a given label (read only). */
  vector<int> &pos(Gene label) const { return (*positions)[label]; }
  /* Reversal operation */
  void reversal(Gene, Gene, IR, IR);
  /* Transposition operation */
  void transposition(Gene, Gene, Gene, IR, IR, IR);
  /* Pretty print genome */
  virtual void serialize(ostream &) const;
};

ostream &operator<<(ostream &os, const Genome &g);
