# MASTER Data Analysis Prototype
This is the data analysis prototype software for MAssively Systematic Transcript End Readout (MASTER) technique. This software analyzes the next generation sequencing results of DNA template library and 5' RNA-Seq in order to obtain the number of RNA reads start at each position of DNA template.

## DNA Template Analysis
The purpose of analyzing the sequencing results of DNA template library is to associate the 7-bp randomized TSS-region with a corresponding second randomized 15-bp barcode region.

C++ program __dna_fastq_parse.cpp__ is used for analyzing sequencing results of DNA template libraries. This program tekes 6 parameters and sequencing results in FASTQ format. The 6 parameters from left to right are sequencing quality score cutoff, length of randomized TSS-region, extra positions need to be considered after randomized TSS-region, length of barcode region, the name of output file containing filtered TSS-regions and barcodes, and the name of output file containing discarded DNA reads. 

The output TSS-regions and barcodes file is tab delimited file with barcode, TSS-region, and count written from left to right. This program will also output a stats file of analyzed DNA template library. 

## 5' RNA-Seq Analysis
The purpose of analyzing 5' RNA-Seq results is to identify the DNA templates and transcription start position of RNA reads.

C++ program __rna_fastq_parse.cpp__ is used for 5' RNA-Seq Analysis. This program takes 8 parameters and sequencing results in FASTQ format. The 8 parameters from left to right are sequencing quality score cutoff, length of digital tag region, length of randomized TSS-region, extra positions need to be considered after randomized TSS-region, length of barcode region, the name of output file containing the summary of RNA reads transcribed from different DNA templates, the name of output stats file of analyzed 5' RNA-Seq results, and the output TSS-regions and barcodes file freom DNA template library analysis. 

Importantly, __rna_fastq_parse.cpp__ will also output a __tag record file__ which records all transcribed RNA TSS-region sequences and their DNA templates and counts. This file is tab delimited with transcribed RNA TSS-region sequences, DNA template TSS-regions, start position at _lacCONS_ template, read counts, number of different digital tags in all reads, and match or mismatch to DNA template TSS-regions.

