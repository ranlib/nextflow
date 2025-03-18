nextflow.enable.dsl=2

include { SampleSV as Sniffles2SV } from './Sniffles2.nf'
include { Clair } from './Clair.nf'

process ZipAndIndexVCF {
    tag "$vcf"

    input:
    path vcf

    output:
    path "${vcf}.gz"
    path "${vcf}.gz.tbi"

    script:
    """
    set -eux
    bgzip -c ${vcf} > ${vcf}.gz
    tabix -p vcf ${vcf}.gz
    """
}

workflow CallVariantsONT {
    params.bam
    params.bai
    params.minsvlen = 50
    params.prefix
    params.sample_id
    params.ref_fasta
    params.ref_fasta_fai
    params.ref_dict
    params.call_svs
    params.call_small_variants
    params.sites_vcf
    params.sites_vcf_tbi

    // Channels for optional outputs
    clair_results = Channel.empty()
    sniffles_results = Channel.empty()

    if (params.call_small_variants) {
        clair_results = Clair(
            params.bam,
            params.bai,
            params.ref_fasta,
            params.ref_fasta_fai,
            params.sites_vcf,
            params.sites_vcf_tbi,
            "ONT"
        )
    }

    if (params.call_svs) {
        sniffles_results = Sniffles2SV(
            params.bam,
            params.bai,
            params.minsvlen,
            params.sample_id,
            params.prefix
        )
        .map { vcf -> ZipAndIndexVCF(vcf) }
    }

    // Define outputs
    output:
    clair_results
    sniffles_results
}
