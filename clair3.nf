nextflow.enable.dsl=2

process Clair {
    tag "clair3_variant_calling"
    
    input:
        path bam
        path bai
        path ref_fasta
        path ref_fasta_fai
        path sites_vcf optional true
        path sites_vcf_tbi optional true
        val chr optional true
        val is_haploid optional true
        val preset
    
    output:
        path "pileup.vcf.gz" optional true
        path "pileup.vcf.gz.tbi" optional true
        path "full_alignment.vcf.gz" optional true
        path "full_alignment.vcf.gz.tbi" optional true
        path "merge_output.vcf.gz" optional true
        path "merge_output.vcf.gz.tbi" optional true
        path "merge_output.gvcf.gz" optional true
        path "merge_output.gvcf.gz.tbi" optional true
    
    script:
        def platform = preset == "CCS" ? "hifi" : "ont"
        """
        set -eux
        SM=$(samtools view -H ${bam} | grep -m1 '^@RG' | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g')
        
        set -euxo pipefail
        
        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)
        
        /opt/bin/run_clair3.sh \${sites_vcf ? "--vcf_fn=${sites_vcf}" : ""} \
        --bam_fn=${bam} \
        --ref_fn=${ref_fasta} \
        --threads=${num_core} \
        --platform=${platform} \
        --model_path="/opt/models/${platform}" \
        --gvcf \
        ${chr ? "--ctg_name=${chr}" : "--include_all_ctgs"} \
        ${is_haploid ? "--no_phasing_for_fa" : ""} \
        --sample_name=${SM} \
        --output="./"
        """
    
    container "hkubal/clair3:v1.0.9"
    cpus 36
    memory "72GB"
    disk "100GB"
}
