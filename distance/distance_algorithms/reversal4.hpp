#pragma once

#include "../misc/dist.hpp"
#include "../misc/permutation.hpp"

class Reversal4 : public DistAlg {
public:
  int estimate_distance(Permutation pi) override;
};
