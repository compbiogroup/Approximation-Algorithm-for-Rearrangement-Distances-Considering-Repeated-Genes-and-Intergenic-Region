#include "ga.hpp"
#include "../misc/io.hpp"
#include <queue>
#include <set>

bool cmp_chr(const unique_ptr<Chromossome> &c1,
             const unique_ptr<Chromossome> &c2) {
  return c1->fitness() > c2->fitness();
}

int GA::select_parent() {
  int index1 = rand() % population_size;
  int index2 = rand() % population_size;
  if ((*population)[index1]->fitness() > (*population)[index2]->fitness()) {
    return index1;
  } else {
    return index2;
  }
}

void GA::select_population(unique_ptr<Population> &mutants) {

  // Find worse mutant.
  vector<int> worse_obj = (*mutants)[0]->fitness();
  int worse_idx = 0;
  for (size_t i = 1; i < mutants->size(); ++i) {
    if ((*mutants)[i]->fitness() < worse_obj) {
      worse_obj = (*mutants)[i]->fitness();
      worse_idx = i;
    }
  }

  // Find best chromosome from the original population.
  vector<int> best_obj = (*population)[0]->fitness();
  int best_idx = 0;
  for (size_t i = 1; i < population->size(); ++i) {
    if ((*population)[i]->fitness() < best_obj) {
      best_obj = (*population)[i]->fitness();
      best_idx = i;
    }
  }

  // Replace worse mutant with best chromosome from original population
  // if the value is better
  if ((*mutants)[worse_idx]->fitness() < (*population)[best_idx]->fitness()) {
    mutants->erase(mutants->begin() + worse_idx);
    mutants->push_back(move((*population)[best_idx]));
  }

  population = move(mutants);
}

bool GA::eval_population(unique_ptr<Population> &new_population, Timer timer) {
  bool improved = false;
  for (auto &chr : *new_population) {
    if (chr->fitness() > best_obj) {
      best_obj = chr->fitness();
      best_chr.reset(new Chromossome(*chr));
      improved = true;
    }
    if (os != nullptr) {
      output(*os, *chr, timer.since_last_mark());
    }
  }
  return improved;
}

/* We include cycles from the original decompositions ignoring conflicts.
 * Afterwards we use bfs to find new cycles. */
Chromossome *GA::crossover(const Chromossome &chr1,
                           const Chromossome &chr2) const {
  auto c_list1 = vector<pair<size_t, Vtx_id>>(chr1.cycle_list());
  auto c_list2 = vector<pair<size_t, Vtx_id>>(chr2.cycle_list());
  Chromossome *chr = new Chromossome(*original);

  random_shuffle(c_list1.begin(), c_list1.end());
  random_shuffle(c_list2.begin(), c_list2.end());
  /* sort(c_list1.begin(), c_list1.end()); */
  /* sort(c_list2.begin(), c_list2.end()); */

  while (c_list1.size() > 0 && c_list2.size() > 0) {
    double p = (double)rand() / RAND_MAX;
    if (p < crossover_rate) {
      chr->check_and_add_cycle(chr1.get_cycle(c_list1.back().second));
      c_list1.pop_back();
    } else {
      chr->check_and_add_cycle(chr2.get_cycle(c_list2.back().second));
      c_list2.pop_back();
    }
  }

  while (c_list1.size() > 0) {
    chr->check_and_add_cycle(chr1.get_cycle(c_list1.back().second));
    c_list1.pop_back();
  }

  while (c_list2.size() > 0) {
    chr->check_and_add_cycle(chr2.get_cycle(c_list2.back().second));
    c_list2.pop_back();
  }

  chr->decompose_with_bfs(true);
  return chr;
}

/* For each cycle we have a chance equals to mutation_rate to removed.
 * Afterwards we use bfs to find new cycles. */
void GA::mutation(unique_ptr<Chromossome> &chr) const {
  for (auto c : chr->cycle_list()) {
    double p = (double)rand() / RAND_MAX;
    if (p < mutation_rate) {
      chr->rem_cycle(c.second);
    }
  }
  chr->decompose_with_bfs(true);
}
