#pragma once

#include "../misc/dist.hpp"
#include "../misc/permutation.hpp"

void transposition_move(Permutation &, Gene, Gene, Gene, IR, IR, IR);

class Transposition3 : public DistAlg {
public:
  int estimate_distance(Permutation pi) override;
};
