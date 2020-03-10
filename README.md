# Modified QuanTest2

A containerized version of QuanTest2.
NOTE: the quantest2.py script is modified to be able to use the standalone version of PSIPRED to predict the secondary structures.

This repository is created to store the Dockerfile and the modified scripts.

The SOURCES are listed below in the REFERENCE section.


# Docker
The docker image is accessible from suzannejin/mod_quantest2
You can also use it with singularity:
singularity pull docker://suzannejin/mod_quantest2:latest

# Usage
```
quantest2 <alignment file> <ss file>
```
The alignment file will be in FASTA format.
The ss file should include the selected 3 reference sequences per family.
You can retrieve the most informative reference sequences by doing:
t_coffee -other_pg seq_reformat -in <ref msa> -action +trim
_aln_n3 -output fasta_seq

# Reference

### Main
QuanTest2:
```
Fabian Sievers, Desmond G Higgins, QuanTest2: benchmarking multiple sequence alignments using secondary structure prediction, Bioinformatics, Volume 36, Issue 1, 1 January 2020, Pages 90–95, https://doi.org/10.1093/bioinformatics/btz552
```
### Dependency
PSIPRED:
```
Jones DT. (1999) Protein secondary structure prediction based on position-specific scoring matrices. J. Mol. Biol. 292: 195-202. 
```
DeepMSA:
```
Chengxin Zhang, Wei Zheng, S M Mortuza, Yang Li, Yang Zhang, DeepMSA: constructing deep multiple sequence alignment to improve contact prediction and fold-recognition for distant-homology proteins, Bioinformatics, , btz863, https://doi.org/10.1093/bioinformatics/btz863
```

[...]






