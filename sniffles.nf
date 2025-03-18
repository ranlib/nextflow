nextflow.enable.dsl=2

process SampleSV {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bai)
    val minsvlen
    val prefix

    output:
    path "${prefix}_${sample_id}.sniffles.snf"
    path "${prefix}_${sample_id}.sniffles.vcf"

    script:
    def snf_output = "${prefix}_${sample_id}.sniffles.snf"
    def vcf_output = "${prefix}_${sample_id}.sniffles.vcf"

    """
    set -eux

    sniffles -t 8 \\
             -i ${bam} \\
             --minsvlen ${minsvlen} \\
             --sample-id ${sample_id} \\
             --vcf ${vcf_output} \\
             --snf ${snf_output}
    """
}

process MergeCall {
    tag "$prefix"

    input:
    path snfs
    val prefix

    output:
    path "${prefix}.vcf"

    script:
    """
    set -eux
    sniffles --input ${snfs.join(' ')} --vcf ${prefix}.vcf
    """
}

workflow Sniffles2 {
    params.sampleBAMs
    params.sampleBAIs
    params.sampleIDs
    params.prefix
    params.minsvlen = 50

    // Run SampleSV for each sample
    sample_results = Channel.fromList(params.sampleIDs)
        .zip(Channel.fromList(params.sampleBAMs), Channel.fromList(params.sampleBAIs))
        .map { sample_id, bam, bai -> tuple(sample_id, bam, bai, params.minsvlen, params.prefix) }
        | SampleSV

    // Collect .snf files and run MergeCall
    sample_results
        .collect()
        | MergeCall(params.prefix)
}
