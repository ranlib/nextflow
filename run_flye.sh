nextflow flye.nf \
--reads ~/Analysis/nanopore/data/S_aureus/reads/S_aureus_JKD6159_ONT_R10.4_guppy_v6.1.7_downsampled_100.fastq.gz \
--prefix testit \
--read_type "\--nano-raw" \
--num_core 2
