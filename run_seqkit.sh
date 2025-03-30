#!/bin/bash
#conda activate analysis

nextflow run seqkit.nf --fastq ERR552797_10k_1.fq.gz --samplename ERR552797 -with-docker "dbest/seqkit:v2.8.2"
