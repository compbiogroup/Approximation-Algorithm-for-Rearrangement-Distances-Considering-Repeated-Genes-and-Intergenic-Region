# Approximation-Algorithm-for-Rearrangement-Distances-Considering-Repeated-Genes-and-Intergenic-Region

Partition and distance algorithms necessary to ensure an Theta(k) approximation for rearrangement distances considering multiple genes and intergenic region sizes as described in the papers "Approximation Algorithm for Rearrangement Distances Considering Repeated Genes and Intergenic Regions" and "Signed Rearrangement Distances Considering Repeated Genes and Intergenic Regions". The `partition` folder contains the partition algorithm and the `distance` folder contains the distance algorithm.

## Instances

Some generated instances are available in the `database` folder. In the file names, `LY` represents the size of the alphabet (`Y`) used to generated the instances.

The following set of instances are present in the folder:
- `REVX_RANDSTR`: unsigned genomes generated with the application of X reversions;
- `SREVX_RANDSTR`: signed genomes generated with the application of X reversions;
- `TRANSX_RANDSTR`: unsigned genomes generated with the application of X transpositions;
- `REVTRANSX_RANDSTR`: unsigned genomes generated with the application of floor of X/2 reversals and ceil of X/2  transpositions;
- `SREVTRANSX_RANDSTR`: signed genomes generated with the application of floor of X/2 reversals and ceil of X/2  transpositions.
