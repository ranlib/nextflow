
process Minimap2 {

    publishDir "alignment", mode: 'copy'
    container "dbest/minimap2:2.28"

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
    """
    set -euxo pipefail
    NUM_CPUS=\$(nproc)
    RAM_IN_GB=\$(free -g | awk '/^Mem:/ {print \$2}')
    MEM_FOR_SORT=\$(echo "" | awk "{print int((\$RAM_IN_GB - 1)/\$NUM_CPUS)}")

    MAP_PARAMS="-aYL --MD --eqx -x ${map_preset} -t \${NUM_CPUS} ${ref_fasta}"
    [[ -n "${RG}" ]] && MAP_PARAMS="\$MAP_PARAMS -R ${RG}"

    SORT_PARAMS="-@\${NUM_CPUS} -m\${MEM_FOR_SORT}G --no-PG -o ${prefix}.pre.bam"
    
    FILE="${reads[0]}"
    FILES="${reads.join(' ')}"

    minimap2 \$MAP_PARAMS ${reads} > tmp.sam
    
    samtools sort \$SORT_PARAMS tmp.sam
    samtools calmd -b --no-PG ${prefix}.pre.bam ${ref_fasta} > ${prefix}.bam
    samtools index -@\${NUM_CPUS} ${prefix}.bam
    """
}