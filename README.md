# Modified QuanTest2

A containerized version of QuanTest2.
NOTE: the _quantest2.py_ script is modified to be able to use the standalone version of PSIPRED to predict the secondary structures.

This repository is created to store the Dockerfile and the modified scripts. A Nextflow pipeline is also included.

The original SOURCES are listed below in the REFERENCE section.


# Docker | Singularity
The docker image is accessible from *suzannejin/mod_quantest2*

You can also use it with singularity:
```
singularity pull docker://suzannejin/mod_quantest2:latest
```

# Usage
```
quantest2 <alignment file> <ss file>
```
The alignment file should be in FASTA format.

The ss file should include the chosen 3 reference sequences per family.
You can retrieve the most informative reference sequences by doing:
```
t_coffee -other_pg seq_reformat -in <ref msa> -action +trim
_aln_n3 -output fasta_seq
```

You can also use the _get_refANDinformative_seqs.py_ script to retrieve N informative sequences from the alignment.

# Nextflow pipeline
The configuration file and the upper part of the pipeline (set input) should be modified before being used, and then run:
```
nextflow run quantest.nf
```

# Reference

## _Main_
__QuanTest2:__

Sievers, F. & Higgins, D. G. QuanTest2: benchmarking multiple sequence alignments using secondary structure prediction. Bioinformatics (2019). doi:10.1093/bioinformatics/btz552

## _Dependency_
__DeepMSA:__

Zhang, C., Zheng, W., Mortuza, S. M., Li, Y. & Zhang, Y. DeepMSA: constructing deep multiple sequence alignment to improve contact prediction and fold-recognition for distant-homology proteins. Bioinformatics (2019). doi:10.1093/bioinformatics/btz863

__PSIPRED:__

Jones, D. T. Protein secondary structure prediction based on position-specific scoring matrices 1 1Edited by G. Von Heijne. J. Mol. Biol. 292, 195–202 (1999).

__T-coffee:__

Notredame, C., Higgins, D. G. & Heringa, J. T-coffee: a novel method for fast and accurate multiple sequence alignment 1 1Edited by J. Thornton. J. Mol. Biol. 302, 205–217 (2000).


[...]






