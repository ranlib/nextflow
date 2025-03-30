#!/bin/bash

READS=~/Analysis/nanopore/data/S_aureus/reads/S_aureus_JKD6159_ONT_R10.4_guppy_v6.1.7_downsampled_100.fastq.gz
FASTA=~/Analysis/nanopore/data/S_aureus/reference/S_aureus_JKD6159.fasta
RG='@RG\\tID:S_aureus_JKD6159\\tSM:S_aureus_JKD6159'

nextflow run alignment.nf --reads $READS --ref_fasta $FASTA --prefix S_aureus --RG $RG --map_preset map-ont