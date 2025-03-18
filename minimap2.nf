nextflow.enable.dsl=2

process Minimap2 {
    tag "$prefix"
    
    input:
    path reads
    path ref_fasta
    val RG
    val map_preset
    val library
    val tags_to_preserve
    val prefix
    
    output:
    path "${prefix}.bam"
    path "${prefix}.bam.bai"
    
    script:
    def num_core = task.cpus ?: 4
    def mem_gb = task.memory ?: 30
    
    def map_params = "-aYL --MD --eqx -x ${map_preset} -t ${num_core} ${ref_fasta}"
    if (RG) {
        map_params += " -R '${RG}'"
    }
    
    """
    set -euxo pipefail
    
    FILES="${reads.join(' ')}"
    if [[ "$FILES" =~ \.(fastq|fq|fasta|fa|fastq.gz|fq.gz|fasta.gz|fa.gz|bam)$ ]]; then
        case "$FILES" in
            *.fastq|*.fq) cat $FILES | minimap2 $map_params - > tmp.sam ;;
            *.fastq.gz|*.fq.gz) zcat $FILES | minimap2 $map_params - > tmp.sam ;;
            *.fasta|*.fa) cat $FILES | python3 /usr/local/bin/cat_as_fastq.py | minimap2 $map_params - > tmp.sam ;;
            *.fasta.gz|*.fa.gz) zcat $FILES | python3 /usr/local/bin/cat_as_fastq.py | minimap2 $map_params - > tmp.sam ;;
            *.bam)
                for f in $FILES; do
                    if [[ -n "${tags_to_preserve}" ]]; then
                        samtools fastq -T ${tags_to_preserve.join(',')} "$f"
                    else
                        samtools fastq "$f"
                    fi
                done > tmp.fastq
                minimap2 $map_params tmp.fastq > tmp.sam ;;
        esac
    else
        echo "Unsupported file format"
        exit 1
    fi
    
    samtools sort -@${num_core} -m${mem_gb}G --no-PG -o ${prefix}.pre.bam tmp.sam
    
    if [[ -n "${library}" ]]; then
        samtools view --no-PG -H ${prefix}.pre.bam > header.txt
        awk '$1 ~ /^@RG/' header.txt > rg_line.txt
        awk -v lib="${library}" '{ for (i=1; i<=NF; ++i) { if ($i ~ "LB:") $i="LB:"lib } print }' rg_line.txt > fixed_rg_line.txt
        sed -n '/@RG/q;p' header.txt > first_half.txt
        sed -n '/@RG/,$p' header.txt | sed '1d' > second_half.txt
        cat first_half.txt fixed_rg_line.txt second_half.txt > fixed_header.txt
        samtools reheader fixed_header.txt ${prefix}.pre.bam > ${prefix}.bam
    else
        samtools calmd -b --no-PG ${prefix}.pre.bam ${ref_fasta} > ${prefix}.bam
    fi
    
    samtools index -@${num_core} ${prefix}.bam
    """
}

workflow Minimap2Workflow {
    params.reads
    params.ref_fasta
    params.RG
    params.map_preset
    params.library = null
    params.tags_to_preserve = []
    params.prefix = "out"
    
    take:
    file reads from params.reads
    file ref_fasta from params.ref_fasta
    val RG from params.RG
    val map_preset from params.map_preset
    val library from params.library
    val tags_to_preserve from params.tags_to_preserve
    val prefix from params.prefix
    
    main:
    Minimap2(reads, ref_fasta, RG, map_preset, library, tags_to_preserve, prefix)
    
    emit:
    bam = Minimap2.out[0]
    bai = Minimap2.out[1]
}
