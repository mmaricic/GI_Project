Illumina paired-end sequencing simulator.

Itroduction
============
Project created as an assignment on course Genome Informatics - School of Electrical Engineering University of Belgrade.
Simulates Illumina paired-end reads created during sequencing.

Instructions for use
=====================
Code is written in jupyter notebook and as such can be used. There is also a Python file that contains the same code.

For running the simulator you need to call "simulatePairedEndReads" function. The parameters are:
- refGenomeFile - path to the FASTA file with reference genome
- quality - average quality of nucleotids. It should be value between 33 and 126
- coverage - Wanted coverage
- readSize - Size of a single read
- insertSize - Size of a fragment
- [OPTIONAL] errorSNV(=0) - Error probability for SNV. Sum of this and errorInDel must not be larger than 1.
- [OPTIONAL] errorInDel(=0) - Error probability for INS/DEL. Sum of this and errorSNV must not be larger than 1.
