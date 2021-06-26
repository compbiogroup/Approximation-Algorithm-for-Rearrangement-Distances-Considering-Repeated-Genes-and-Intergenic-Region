# Approximation for Minimum Common Intergenic String Partition

An $2k$ approximation algorithm for the Minimum Common Intergenic String Partition and Reverse Minimum Common Intergenic String Partition, where $k$ in the maximum number of copies of a given character in the strings.

## Usage

Compile the code by running `make` and see the running options with `./MCISP2K --help`. A program to produce test instances is also available, see `./DB --help`.

## Instances

Some generated instances are available in `database`. Each file in `REVX_RANDSTR` was generated with the aplication of X reversions, each file in `TRANSX_RANDSTR` was generated with the application of X transpositions, and each file in `REVTRANSX_RANDSTR` was generated with the application of floor of X/2 reversals and ceil of X/2  transpositions. In the file names, `LX` represents the size of the alphabet used to generated the instances.
