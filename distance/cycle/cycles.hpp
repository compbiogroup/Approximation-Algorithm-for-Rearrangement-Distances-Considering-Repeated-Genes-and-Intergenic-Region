#pragma once

#include "../misc/genome.hpp"
#include <bitset>
#include <map>
#include <queue>
#include <set>

typedef int Vtx_id;

struct Vertex {
  Vtx_id black;
  int weigth;
  Vtx_id fix_gray;
  vector<Vtx_id> grays;
  bool in_cycle;
  bool sign_positive;
  Vertex() : grays() {
    black = -1;
    weigth = 0;
    fix_gray = -1;
    in_cycle = false;
  }
};

struct PermsIrs {
  size_t n;
  int *s;
  int *s_ir;
  int *p;
  int *p_ir;
};

struct QEntry {
  Vtx_id vtx;
  map<Vtx_id, Vtx_id> fixed;
  set<Vtx_id> vizited;

  QEntry(Vtx_id start) : fixed(), vizited() { vtx = start; }

  QEntry(Vtx_id v, Vtx_id u, Vtx_id alter_v, Vtx_id alter_u,
         map<Vtx_id, Vtx_id> fixed_old, set<Vtx_id> vizited_old)
      : fixed(fixed_old), vizited(vizited_old) {
    vtx = u;
    fixed[v] = u;
    fixed[u] = v;
    fixed[alter_v] = alter_u;
    fixed[alter_u] = alter_v;
    vizited.insert(u);
  }

  QEntry(Vtx_id v, Vtx_id u, map<Vtx_id, Vtx_id> fixed_old,
         set<Vtx_id> vizited_old)
      : fixed(fixed_old), vizited(vizited_old) {
    vtx = u;
    fixed[v] = u;
    fixed[u] = v;
    vizited.insert(u);
  }
};

class CycleGraph {

private:
  int balanced_cycles = 0;
  vector<Vertex> vertices;
  vector<pair<size_t, Vtx_id>>
      cycles; // we indentify cycles by one of their vertices and their sizes
  void bfs(Vtx_id start, unique_ptr<queue<QEntry>> &q,
           unique_ptr<set<pair<Vtx_id, int>>> &vizited_in_level,
           bool is_random);

public:
  /* Initial Constructor. */
  CycleGraph(const Genome &origin, const Genome &target);
  /* Copy Constructor. */
  CycleGraph(const CycleGraph &that)
      : vertices(that.vertices), cycles(that.cycles) {}
  size_t size() const { return vertices.size(); };
  void decompose_with_bfs(bool is_random);
  /* Select a cycle with a bfs
   * Arguments:
   *     start - initial vertex
   *     is_random - whether to use random approach
   */
  /* Get string with the cycles from the decomposition */
  void bfs(Vtx_id start, bool is_random);
  string show_cycles() const;
  /* Recover cycles from string */
  void read_cycles(string);
  /* Recover permutations from decomposition */
  PermsIrs get_perms();
  /* Get number of cycles from the decomposition */
  int dec_size() const { return cycles.size(); }
  int dec_balanced_cycles() const {return balanced_cycles; }
  void rem_cycle(Vtx_id i);
  /* Verify if a cycle can be added */
  bool check_cycle(vector<Vtx_id> cycle);
  void add_cycle(vector<Vtx_id> cycle);
  void check_and_add_cycle(vector<Vtx_id> cycle);
  vector<pair<size_t, Vtx_id>> cycle_list() const { return cycles; }
  vector<Vtx_id> get_cycle(Vtx_id i) const;
  int cycle_weight(Vtx_id i) const;
  void serialize(ostream &) const;
};

ostream &operator<<(std::ostream &os, const CycleGraph &cg);
