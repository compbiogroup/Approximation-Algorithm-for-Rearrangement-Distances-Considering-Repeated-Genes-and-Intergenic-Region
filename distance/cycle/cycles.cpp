#include "cycles.hpp"
#include "../misc/genome.hpp"
#include "../misc/io.hpp"
#include "greedy-dup-cycles.hpp"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <numeric>
#include <queue>
#include <sstream>

/* Helper function to find the edges */
void fill_edges(Genome &s1, Genome &s2, vector<Vertex> &vec) {}

CycleGraph::CycleGraph(const Genome &origin, const Genome &target)
    : vertices(2 * (origin.size() + target.size()) - 4), cycles() {
  /* Black Edges */
  for (size_t i = 0; i < vertices.size(); ++i) {
    if (i % 2 == 0) {
      vertices[i].black = i + 1;
    } else {
      vertices[i].black = i - 1;
    }
  }

  int fhs = vertices.size() / 2; // size of first half of the vertices;

  /* Weigths */
  vertices[0].weigth = origin.get_ir(1);
  for (size_t i = 1; i < origin.size() - 1; ++i) {
    vertices[2 * i - 1].weigth = origin.get_ir(i);
    vertices[2 * i].weigth = origin.get_ir(i + 1);
  }
  vertices[fhs - 1].weigth = origin.get_ir(origin.size() - 1);
  vertices[fhs].weigth = target.get_ir(1);
  for (size_t i = 1; i < target.size() - 1; ++i) {
    vertices[fhs + 2 * i - 1].weigth = target.get_ir(i);
    vertices[fhs + 2 * i].weigth = target.get_ir(i + 1);
  }
  vertices[vertices.size() - 1].weigth = target.get_ir(target.size() - 1);

  /* Signed */
  vertices[0].sign_positive = true;
  for (size_t i = 1; i < origin.size() - 1; ++i) {
    vertices[2 * i - 1].sign_positive = origin[i + 1] >= 0;
    vertices[2 * i].sign_positive = origin[i + 1] >= 0;
  }
  vertices[fhs - 1].sign_positive = true;
  vertices[fhs].sign_positive = true;
  for (size_t i = 1; i < target.size() - 1; ++i) {
    vertices[fhs + 2 * i - 1].sign_positive = target[i + 1] >= 0;
    vertices[fhs + 2 * i].sign_positive = target[i + 1] >= 0;
  }
  vertices[vertices.size() - 1].sign_positive = true;

  /* Gray Edges */
  assert(origin[1] == target[1] &&
         origin[origin.size()] == target[target.size()] &&
         origin.pos(origin[1]).size() == 1 &&
         origin.pos(origin[origin.size()]).size() == 1 &&
         target.pos(target[1]).size() == 1 &&
         target.pos(target[target.size()]).size() == 1);

  vertices[0].grays.push_back(fhs);
  vertices[fhs].grays.push_back(0);
  vertices[fhs - 1].grays.push_back(vertices.size() - 1);
  vertices[vertices.size() - 1].grays.push_back(fhs - 1);
  for (size_t i = 2; i < origin.size(); ++i) {
    for (auto j : target.pos(abs(origin[i]))) {
      if (origin[i] == target[j]) {
        vertices[2 * i - 3].grays.push_back(fhs + 2 * j - 3);
        vertices[2 * i - 2].grays.push_back(fhs + 2 * j - 2);
        vertices[fhs + 2 * j - 3].grays.push_back(2 * i - 3);
        vertices[fhs + 2 * j - 2].grays.push_back(2 * i - 2);
      } else {
        vertices[2 * i - 2].grays.push_back(fhs + 2 * j - 3);
        vertices[2 * i - 3].grays.push_back(fhs + 2 * j - 2);
        vertices[fhs + 2 * j - 2].grays.push_back(2 * i - 3);
        vertices[fhs + 2 * j - 3].grays.push_back(2 * i - 2);
      }
    }
  }
  for (size_t i = 0; i < vertices.size(); ++i) {
    if (vertices[i].grays.size() == 1) {
      vertices[i].fix_gray = vertices[i].grays[0];
    }
  }
}

/* Decompose the remaning graph using bfs */
void CycleGraph::decompose_with_bfs(bool is_random) {
  vector<Vtx_id> idxs(size());
  iota(idxs.begin(), idxs.end(), 0);
  if (is_random) {
    random_shuffle(idxs.begin(), idxs.end());
  }
  for (Vtx_id i : idxs) {
    bfs(i, is_random);
  }
}

void CycleGraph::bfs(Vtx_id start, bool is_random) {
  unique_ptr<QEntry> entry;
  vector<unique_ptr<QEntry>> q1;
  vector<unique_ptr<QEntry>> q2;
  set<Vtx_id> vizited_in_level;
  Vtx_id headtailcorresp_v, headtailcorresp_u;
  size_t pos = 0;

  if (vertices[start].in_cycle)
    return;

  q1.push_back(unique_ptr<QEntry>(new QEntry(start)));
  vizited_in_level.insert(start);

  do {
    entry = move(q1[pos]);
    pos++;

    Vtx_id v = vertices[entry->vtx].black;
    vector<Vtx_id> neigs;
    if (vertices[v].fix_gray != -1) {
      neigs.push_back(vertices[v].fix_gray);
    } else {
      auto vu = entry->fixed.find(v);
      if (vu != entry->fixed.end()) {
        neigs.push_back(vu->second);
      } else {
        for (auto u : vertices[v].grays) {
          if (entry->fixed.find(u) == entry->fixed.end() &&
              vertices[u].fix_gray == -1) {
            neigs.push_back(u);
          }
        }
      }
    }

    vector<Vtx_id> notDone;
    for (Vtx_id u : neigs) {
      if (vizited_in_level.find(u) == vizited_in_level.end() &&
          entry->vizited.find(u) == entry->vizited.end()) {
        notDone.push_back(u);
      }
    }

    if (vertices[v].fix_gray == -1) {
      for (Vtx_id u : notDone) {
        if ((v % (vertices.size() / 2)) % 2 == 1) {
          headtailcorresp_v = v + 1;
        } else {
          headtailcorresp_v = v - 1;
        }
        if ((u % (vertices.size() / 2)) % 2 == 1) {
          headtailcorresp_u = u + 1;
        } else {
          headtailcorresp_u = u - 1;
        }
        q2.push_back(unique_ptr<QEntry>(
            new QEntry(v, u, headtailcorresp_v, headtailcorresp_u, entry->fixed,
                       entry->vizited)));
        vizited_in_level.insert(u);
      }
    } else {
      for (Vtx_id u : notDone) {
        q2.push_back(
            unique_ptr<QEntry>(new QEntry(v, u, entry->fixed, entry->vizited)));
        vizited_in_level.insert(u);
      }
    }

    if (pos == q1.size()) {
      q1.clear();
      for (unique_ptr<QEntry> &e : q2) {
        q1.push_back(move(e));
      }
      if (is_random) {
        random_shuffle(q1.begin(), q1.end());
      }
      q2.clear();
      vizited_in_level.clear();
      pos = 0;
    }
  } while (q1[pos]->vtx != start);

  /* Find balanced cycle in level if it exists */
  vector<Vtx_id> cycle;
  bool ok = false;
  for (;pos < q1.size(); pos++) {
      if (q1[pos]->vtx == start) {
          entry = move(q1[pos]);
          cycle.clear();
          Vtx_id u, v;
          int weigth = 0;
          v = start;
          do {
            weigth += vertices[v].weigth;
            u = vertices[v].black;
            cycle.push_back(v);
            cycle.push_back(u);
            v = entry->fixed[u];
          } while (v != start);
          if (weigth == 0) {
              ok = true;
              add_cycle(cycle);
              break;
          }
      }
  }
  if (!ok) {
      add_cycle(cycle);
  }
}

string CycleGraph::show_cycles() const {
  ostringstream ss;

  ss << '[';
  for (size_t i = 0; i < cycles.size(); ++i) {
    Vtx_id v = cycles[i].second;
    ss << '[' << v;
    v = vertices[v].black;
    ss << "," << v;
    v = vertices[v].fix_gray;
    while (v != cycles[i].second) {
      ss << "," << v;
      v = vertices[v].black;
      ss << "," << v;
      v = vertices[v].fix_gray;
    }
    ss << ']';
    if (i != cycles.size() - 1)
      ss << ',';
  }
  ss << ']';

  return ss.str();
}

void CycleGraph::read_cycles(string str) {
  string str_cycle, str_el;
  vector<Vtx_id> cycle;

  /* Remove '[[' and ']]'. */
  str = str.substr(2, str.size() - 4);

  /* Read each cycle. */
  string delimiter = "],[";
  while (not str.empty()) {
    size_t idx = str.find(delimiter);
    if (idx != str.npos) {
      str_cycle = str.substr(0, idx);
      str.erase(0, idx + delimiter.length());
    } else {
      str_cycle = str;
      str = "";
    }

    /* Read each element. */
    cycle.clear();
    stringstream ss(str_cycle);
    while (getline(ss, str_el, ',')) {
      cycle.push_back(stoi(str_el));
    }

    add_cycle(cycle);
  }
}

PermsIrs CycleGraph::get_perms() {
  size_t n = (size() - 4) / 4 + 2;
  PermsIrs perms_irs;
  perms_irs.n = n;
  perms_irs.s = (int *)malloc(n * sizeof(int));
  perms_irs.s_ir = (int *)malloc((n - 1) * sizeof(int));
  perms_irs.p = (int *)malloc(n * sizeof(int));
  perms_irs.p_ir = (int *)malloc((n - 1) * sizeof(int));

  /* The second permutation is the identity. */
  for (size_t i = 0; i < n; ++i) {
    perms_irs.p[i] = i;
  }

  for (size_t i = 0; i < cycles.size(); ++i) {
    Vtx_id v = cycles[i].second;
    Vtx_id u = vertices[v].black;
    v = vertices[u].fix_gray;
    Vtx_id a = min(v, u);
    Vtx_id b = max(v, u);
    if (vertices[a].sign_positive == vertices[b].sign_positive) {
      perms_irs.s[(a + 1) / 2] = perms_irs.p[(b - size() / 2 + 1) / 2];
    } else {
      perms_irs.s[(a + 1) / 2] = -perms_irs.p[(b - size() / 2 + 1) / 2];
    }
    while (v != cycles[i].second) {
      u = vertices[v].black;
      v = vertices[u].fix_gray;
      Vtx_id a = min(v, u);
      Vtx_id b = max(v, u);
      if (vertices[a].sign_positive == vertices[b].sign_positive) {
        perms_irs.s[(a + 1) / 2] =
            perms_irs.p[(b - size() / 2 + 1) / 2];
      } else {
        perms_irs.s[(a + 1) / 2] =
            -perms_irs.p[(b - size() / 2 + 1) / 2];
      }
    }
  }

  for (size_t i = 1; i < n; ++i) {
    perms_irs.s_ir[i-1] = vertices[2*i - 1].weigth;
    perms_irs.p_ir[i-1] = vertices[vertices.size() / 2 + 2*i - 1].weigth;
  }

  return perms_irs;
}

void CycleGraph::serialize(ostream &os) const {
  for (size_t i = 0; i < vertices.size(); ++i) {
    os << i << "(" << vertices[i].black << "," << vertices[i].weigth << "," << vertices[i].in_cycle << "," << vertices[i].sign_positive << "|";
    os << "[";
    if (vertices[i].fix_gray == -1) {
      for (auto gray : vertices[i].grays) {
        os << gray << ",";
      }
    } else {
      os << vertices[i].fix_gray;
    }
    os << "])"
       << " ";
    if (i == (vertices.size() - 1) / 2)
      os << endl;
  }
}

ostream &operator<<(std::ostream &os, const CycleGraph &cg) {
  cg.serialize(os);
  return os;
}

void CycleGraph::rem_cycle(Vtx_id i) {
  Vtx_id headtailcorresp_v, headtailcorresp_u;

  if (cycle_weight(i) == 0) {
      balanced_cycles--;
  }
  vector<Vtx_id> cycle = get_cycle(i);

  assert(cycle.size() % 2 == 0);
  for (size_t i = 0; i < cycle.size(); ++i) {
    Vtx_id v = cycle[i];
    assert(vertices[v].in_cycle);
    vertices[v].in_cycle = false;
    if (i % 2 == 1 && vertices[v].grays.size() > 1) {
      Vtx_id u = (i == cycle.size() - 1) ? cycle[0] : cycle[i + 1];
      if ((v % (vertices.size() / 2)) % 2 == 1) {
        headtailcorresp_v = v + 1;
      } else {
        headtailcorresp_v = v - 1;
      }
      if ((u % (vertices.size() / 2)) % 2 == 1) {
        headtailcorresp_u = u + 1;
      } else {
        headtailcorresp_u = u - 1;
      }
      if (not vertices[headtailcorresp_v].in_cycle && not vertices[headtailcorresp_u].in_cycle) {
        assert(vertices[v].fix_gray != -1);
        assert(vertices[u].fix_gray != -1);
        assert(vertices[headtailcorresp_v].fix_gray != -1);
        assert(vertices[headtailcorresp_u].fix_gray != -1);
        vertices[v].fix_gray = -1;
        vertices[u].fix_gray = -1;
        vertices[headtailcorresp_v].fix_gray = -1;
        vertices[headtailcorresp_u].fix_gray = -1;
      }
    }
  }
  cycles.erase(find(cycles.begin(), cycles.end(),
                    pair<size_t, Vtx_id>(cycle.size(), cycle[0])));
}

bool CycleGraph::check_cycle(vector<Vtx_id> cycle) {
  for (size_t i = 0; i < cycle.size(); ++i) {
    Vtx_id v = cycle[i];
    if (vertices[v].in_cycle)
      return false;
    if (i % 2 == 1) {
      Vtx_id u = (i == cycle.size() - 1) ? cycle[0] : cycle[i + 1];
      if (vertices[v].fix_gray != -1) {
        if (vertices[v].fix_gray != u)
          return false;
      } else if (vertices[u].fix_gray != -1) {
        return false;
      }
    }
  }
  return true;
}

void CycleGraph::add_cycle(vector<Vtx_id> cycle) {
  Vtx_id headtailcorresp_v, headtailcorresp_u;

  assert(cycle.size() % 2 == 0);
  for (size_t i = 0; i < cycle.size(); ++i) {
    Vtx_id v = cycle[i];
    assert(not vertices[v].in_cycle);
    vertices[v].in_cycle = true;
    if (i % 2 == 1 && vertices[v].fix_gray == -1) {
      Vtx_id u = (i == cycle.size() - 1) ? cycle[0] : cycle[i + 1];
      assert(vertices[u].fix_gray == -1);
      vertices[v].fix_gray = u;
      vertices[u].fix_gray = v;
      if ((v % (vertices.size() / 2)) % 2 == 1) {
        headtailcorresp_v = v + 1;
      } else {
        headtailcorresp_v = v - 1;
      }
      if ((u % (vertices.size() / 2)) % 2 == 1) {
        headtailcorresp_u = u + 1;
      } else {
        headtailcorresp_u = u - 1;
      }
      assert(vertices[headtailcorresp_u].fix_gray == -1);
      assert(vertices[headtailcorresp_v].fix_gray == -1);
      vertices[headtailcorresp_v].fix_gray = headtailcorresp_u;
      vertices[headtailcorresp_u].fix_gray = headtailcorresp_v;
    }
  }
  cycles.push_back(pair<size_t, Vtx_id>(cycle.size(), cycle[0]));
  if (cycle_weight(cycle[0]) == 0) {
      balanced_cycles++;
  }
}

void CycleGraph::check_and_add_cycle(vector<Vtx_id> cycle) {
  if (check_cycle(cycle)) {
    add_cycle(cycle);
  }
}

vector<Vtx_id> CycleGraph::get_cycle(Vtx_id i) const {
  vector<Vtx_id> cycle;

  Vtx_id v = i;
  cycle.push_back(v);
  v = vertices[v].black;
  cycle.push_back(v);
  v = vertices[v].fix_gray;
  while (v != cycle[0]) {
    cycle.push_back(v);
    v = vertices[v].black;
    cycle.push_back(v);
    v = vertices[v].fix_gray;
  }

  return cycle;
}

int CycleGraph::cycle_weight(Vtx_id i) const {
  vector<Vtx_id> cycle;
  int weigth = 0;

  Vtx_id v = i;
  weigth += vertices[v].weigth;
  v = vertices[v].black;
  v = vertices[v].fix_gray;
  while (v != i) {
    weigth += vertices[v].weigth;
    v = vertices[v].black;
    v = vertices[v].fix_gray;
  }

  return weigth;
}
