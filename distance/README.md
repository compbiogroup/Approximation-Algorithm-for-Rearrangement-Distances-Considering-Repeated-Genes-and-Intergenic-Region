# 3 Asymptotic Approximation for Sorting by Intergenic Transposition

A 3 asymptotic approximation algorithm for Sorting by Intergenic Transposition. The algorithm can be applied to genomes with multiple genes, in that case random mappings of the genomes into permutations will be used to produce the distances.

## Usage

Compile the code by running `make` and see the running options with `./dist --help`.  To use other algorithm besides the transposition approximation include a executable in the `external` folder and pass its name as the algorithm. In that case the executable should receive one instance via the standard input and produce the distance in the standard output.
