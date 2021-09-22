#pragma once

#include "../misc/dist.hpp"
#include "../misc/permutation.hpp"

class ReversalTransposition45 : public DistAlg {
public:
  int estimate_distance(Permutation pi) override;
};
