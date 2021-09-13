#include "Transposition4/transposition4.hpp"
#include "external/external.hpp"
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
  string input_file;
  string output_folder;
  int iterations = 1;
  bool extend = false;
  bool duplicate = false;
  string alg;
};

void help(char *name) {
  cout << "usage: Calculate distances for pairs of strings. If replicas are "
       << "presents multiple random mappings are generated." << endl
       << "\t" << name << " ALG [OPTIONS]" << endl
       << endl
       << "positional arguments:" << endl
       << "\tALG                     the algorithm to use (transposition4, or "
          "the name of an executable in the external folder)"
       << endl
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
       << "\t-k, --iterations ITER   number of iterations (defaut 1)" << endl
       << endl
       << "\t-e, --extend            whether to extend the genomes before "
          "apply the algorithm"
       << endl
       << "\t-d, --duplicate         whether to duplicate sign genes to turn "
          "then "
          "into unsigned genes"
       << endl;

  exit(EXIT_SUCCESS);
}

void get_args(Args &args, int argc, char *argv[]) {
  extern char *optarg;
  extern int optind;
  int n_pos_args = 0;

  struct option longopts[] = {
      {"input", 1, NULL, 'i'},      {"output", 1, NULL, 'o'},
      {"iterations", 1, NULL, 'k'}, {"extend", 0, NULL, 'e'},
      {"duplicate", 0, NULL, 'd'},  {"help", 0, NULL, 'h'}};

  char op;
  while ((op = getopt_long(argc, argv, "i:o:k:hed", longopts, NULL)) != -1) {
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
    case 'd':
      args.duplicate = true;
      break;
    default:
      help(argv[0]);
    }
  }
  for (int i = optind; i < argc; i++) {
    args.alg = argv[i];
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

  // set algorithm
  if (args.alg == "transposition4") {
    alg.reset(new Transposition4());
  } else {
    alg.reset(new ExternalDistAlg("external/" + args.alg));
  }

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
      int dist, dist_best = std::numeric_limits<int>::max();

      InputData data =
          input((*input_lines)[i], (*input_lines)[i + 1], (*input_lines)[i + 2],
                (*input_lines)[i + 3], args.extend);
      auto g = unique_ptr<Genome>(data.g);
      auto h = unique_ptr<Genome>(data.h);

      if (args.output_folder != "") {
        os.open(args.output_folder + '/' +
                fs::path(args.input_file).filename().c_str() +
                string(5 - to_string(i / 4).size(), '0') + to_string(i / 4) +
                "-all");
      }

      for (int j = 1; j <= args.iterations; ++j) {
        timer.mark_time();
        pi.reset(new Permutation(*g, *h, args.duplicate));
        dist = alg->estimate_distance(*pi);
        output((args.output_folder != "") ? os : cout, dist,
               timer.since_last_mark());
        if (dist < dist_best) {
          dist_best = dist;
          pi_best = move(pi);
        }
      }

      if (args.output_folder != "") {
        os.close();
        os.open(args.output_folder + '/' +
                fs::path(args.input_file).filename().c_str() +
                string(5 - to_string(i / 4).size(), '0') + to_string(i / 4) +
                "-best");
      }

      if (args.output_folder != "") {
        os << *pi_best << endl;
        output(os, dist_best, timer.elapsed_time());
      } else {
        cout << *pi_best << endl;
        output(cout, dist_best, timer.elapsed_time());
      }
    }

  } catch (const invalid_argument &e) {
    cerr << "Something went wrong!!!" << endl;
    cerr << e.what() << endl;
  }

  return 0;
}
