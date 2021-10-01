# Approximation Algorithms to Rearrangement Distance on Permutations

The following algorithms are implemented:
- A 3 asymptotic approximation algorithm for Sorting by Intergenic Transposition.
- An 4-approximation algorithm for Sorting by Intergenic Reversals [[1]](#1).
- An 4.5-approximation algorithm for Sorting by Intergenic Reversals and Transposition [[1]](#1).

The algorithms can be applied to genomes with multiple genes, in that case random mappings of the genomes into permutations will be used to produce the distances.

## Usage

Compile the code by running `make` and see the running options with `./dist --help`.  To use other algorithm when calculating the distances include a executable in the `external` folder and pass its name as the algorithm parameter. In that case, the executable should receive one instance as command line arguments (four coma separated lists) and produce the distance in the standard output.

## References

<a id="1">[1]</a> 
Klairton Lima Brito, Géraldine Jean, Guillaume Fertin, Andre Rodrigues Oliveira, Ulisses Dias and Zanoni Dias. Sorting by genome rearrangements on both gene order and intergenic sizes. Journal of Computational Biology 2020;27(2):156–74.
