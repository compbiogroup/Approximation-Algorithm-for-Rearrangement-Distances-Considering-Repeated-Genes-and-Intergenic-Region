#include "cycle/cycles.hpp"
#include "external/external.hpp"
#include "heur/ga.hpp"
#include "misc/genome.hpp"
#include "misc/io.hpp"
#include "misc/permutation.hpp"
#include "misc/timer.hpp"
#include <experimental/filesystem>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
namespace fs = experimental::filesystem;
using namespace std;

#define N_POS_ARGS 1

struct Args {
  string heuristic;
  string input_file;
  string output_folder;
  int iterations = 1;
  bool extend = false;
  bool duplicate = false;
  bool all = false;
};

void help(char *name) {
  cout << "usage: Find a cicle decomposition for the breakpoint graph formed "
          "be the origin and target signed genome."
       << endl
       << "\t" << name << " HEUR [OPTIONS]" << endl
       << endl
       << "positional arguments:" << endl
       << "\tHEUR                    the heuristic to use (ga|rand)" << endl
       << endl
       << "optional arguments:" << endl
       << "\t-h, --help              show this help message and exit" << endl
       << "\t-i, --input INPUT       input file (if not provided stdin is "
          "used). Each 4 lines of the input file correspond to a instance, "
          "each line has a list of space separated values, and represent in "
          "order the origin string, the origin intergenic region list, the "
          "target string, and the target intergenic region list."
       << endl
       << "\t-o, --output OUTPUT     outpute folder (if not provided stdout is "
          "used)"
       << endl
       << "\t-k, --iterations ITER   number of iterations (default 1)" << endl
       << endl
       << "\t-e, --extend            whether to extend the genomes before "
          "apply the algorithm"
       << endl
       << "\t-a, --all               print all cycles" << endl
       << endl;

  exit(EXIT_SUCCESS);
}

void get_args(Args &args, int argc, char *argv[]) {
  extern char *optarg;
  extern int optind;
  int n_pos_args = 0;

  struct option longopts[] = {{"input", 1, NULL, 'i'},
                              {"output", 1, NULL, 'o'},
                              {"iterations", 1, NULL, 'k'},
                              {"extend", 0, NULL, 'e'},
                              {"help", 0, NULL, 'h'}};

  char op;
  while ((op = getopt_long(argc, argv, "i:o:k:heda", longopts, NULL)) != -1) {
    switch (op) {
    case 'i':
      args.input_file = optarg;
      break;
    case 'o':
      args.output_folder = optarg;
      break;
    case 'k':
      args.iterations = atoi(optarg);
      break;
    case 'e':
      args.extend = true;
      break;
    case 'a':
      args.all = true;
      break;
    default:
      help(argv[0]);
    }
  }
  for (int i = optind; i < argc; i++) {
    args.heuristic = argv[i];
    n_pos_args++;
  }

  if (n_pos_args != N_POS_ARGS) {
    help(argv[0]);
  }
}

int main(int argc, char *argv[]) {
  Args args;
  ifstream is;
  unique_ptr<vector<string>> input_lines;
  unique_ptr<DistAlg> alg;

  get_args(args, argc, argv);

  // set seed
  srand(1);

  if (args.input_file != "") {
    is.open(args.input_file);
    input_lines.reset(read_lines(is));
    is.close();
  } else {
    input_lines.reset(read_lines(cin));
  }

  try {
    if (input_lines->size() % 4 == 1) {
      throw invalid_argument("Number of lines is not multiple of 4.");
    }
#pragma omp parallel for
    for (size_t i = 0; i < input_lines->size(); i += 4) {
      Timer timer;
      ofstream os;
      unique_ptr<Permutation> pi, pi_best;
      unique_ptr<CycleGraph> cg, cg_aux, cg_best;
      int name_idx = i / 4;

      InputData data =
          input((*input_lines)[i], (*input_lines)[i + 1], (*input_lines)[i + 2],
                (*input_lines)[i + 3], args.extend);
      cg = unique_ptr<CycleGraph>(new CycleGraph(*data.g, *data.h));

      if (args.output_folder != "" && args.all) {
        os.open(args.output_folder + '/' +
                fs::path(args.input_file).filename().c_str() +
                string(5 - to_string(name_idx).size(), '0') +
                to_string(name_idx) + "-all");
      }

      if (args.heuristic == "rand") {
        cg_aux.reset(new CycleGraph(*cg));
        cg_aux->decompose_with_bfs(false);
        cg_best = move(cg_aux);
        if (args.all) {
          output((args.output_folder != "") ? os : cout, *cg_aux,
                 timer.since_last_mark());
        }
        for (int i = 1; i < args.iterations; ++i) {
          timer.mark_time();
          cg_aux.reset(new CycleGraph(*cg));
          cg_aux->decompose_with_bfs(true);
          if (args.all) {
            output((args.output_folder != "") ? os : cout, *cg_aux,
                   timer.since_last_mark());
          }
          if (cg_aux->dec_size() > cg_best->dec_size()) {
            cg_best = move(cg_aux);
          }
        }
      } else if (args.heuristic == "ga") {
        ostream *ga_os;
        if (args.all) {
          ga_os = (args.output_folder != "") ? &os : &cout;
        } else {
          ga_os = nullptr;
        }
        int start = args.iterations / 10;
        if (start % 2 == 1) {
          start += 1;
        }
        GA ga = GA(new Chromossome(*cg), 0.5, 0.5, start, start,
                   (args.iterations - start) / start, ga_os, timer);
        ga.solve(timer);
        cg_best = ga.get_best_chr();
      } else {
        help(argv[0]);
      }

      if (args.output_folder != "") {
        os.close();
        os.open(args.output_folder + '/' +
                fs::path(args.input_file).filename().c_str() +
                string(5 - to_string(name_idx).size(), '0') +
                to_string(name_idx) + "-best");
      }

      output((args.output_folder != "") ? os : cout, *cg_best,
             timer.elapsed_time());

      if (args.output_folder != "") {
        os.close();
        os.open(args.output_folder + '/' +
                fs::path(args.input_file).filename().c_str() +
                string(5 - to_string(name_idx).size(), '0') +
                to_string(name_idx) + "-perm");
      }

      output((args.output_folder != "") ? os : cout, cg_best->get_perms());
    }

  } catch (const invalid_argument &e) {
    cerr << "Something went wrong!!!" << endl;
    cerr << e.what() << endl;
  }

  return 0;
}
