#pragma once
#include "genome.hpp"
using namespace std;

class Permutation : public Genome {
  unique_ptr<vector<IR>> intergenic_regions_target;

public:
  Permutation() {};
  /* Copy constructor */
  Permutation(const Permutation &pi);
  /* Construct a permutation replicas are mapped randomly */
  Permutation(const Genome &s, const Genome &h, bool duplicate);
  /* Verify if the permutation is sorted and with the correct intergenic region
   * sizes */
  bool is_iota() const;
  /* Only read access to intergenic regions of target. */
  IR get_ir_target(int i) const { return (*intergenic_regions_target)[i - 1]; }
  bool breakpoint(Gene) const;
  bool hard_breakpoint(Gene) const;
  bool soft_breakpoint(Gene) const;
  bool overcharged_breakpoint(Gene) const;
  bool undercharged_breakpoint(Gene) const;
  virtual void serialize(ostream &) const override;
};
