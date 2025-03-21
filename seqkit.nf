process stats {
    publishDir "seqkit", mode: 'copy'
    container "dbest/seqkit:v2.8.2"

    input:
    path fastq
    val samplename 

    output:
    path "${samplename}.stats"

    script:
    """
    seqkit stats -T -a -o "${samplename}.stats" $fastq
    """
}

workflow {
    stats_out = stats(file(params.fastq), params.samplename)
    stats_out.view{"${it.text}"}
}
