Illumina paired-end read sequencing simulator.

Itroduction
============
Project created as an assignment on course Genome Informatics - School of Electrical Engineering University of Belgrade.
Simulates Illumina paired-end reads created during sequencing.

Content
==============
**Simulator** - folder containing Simulator and Tests for it written in Jupyter notebook
  SequencingSimulator.ipynb - Simulator for paired-end reads
  Tests.ipynb - Tests for simulator
**BWA-MEM Evaluation** - folder contaning neseccary python scripts for testing BWA-MEM tool
  SequencingSimulator.py - Simulator for paired-end reads (clean python)
  BWA-MEM_AlignmentEvaluation.py - Script for evaluating quality of alignments for BWA-MEM tool based on SAM file produced by our simulator
**Examples** - folder containing examples for testing the alghorithm

Instructions for use
=====================
Code is written in jupyter notebook and as such can be used.

For running the simulator you need to call "simulatePairedEndReads" function. The parameters are:
- refGenomeFile - path to the FASTA file with reference genome
- quality - average quality of nucleotids. It should be value between 33 and 126
- coverage - Wanted coverage
- readSize - Size of a single read
- insertSize - Size of a fragment
- [OPTIONAL] errorSNV(=0) - Error probability for SNV. Sum of this and errorInDel must not be larger than 1.
- [OPTIONAL] errorInDel(=0) - Error probability for INS/DEL. Sum of this and errorSNV must not be larger than 1.

As a result, 2 FASTQ files containg paired-end reads are created. Also it generates SAM file containing the alignments of reads in reference genome.
