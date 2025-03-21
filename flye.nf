process Assemble {
    tag "$prefix"
    
    publishDir "flye", mode: 'copy'
    container "dbest/flye:v2.9.5"

    input:
    path reads
    val prefix
    val read_type
    
    output:
    path "${prefix}.flye.fa"
    path "${prefix}.flye.gfa"
    
    script:
    def num_core = System.getenv('NXF_EXECUTOR') == 'local' ? Runtime.runtime.availableProcessors() : task.cpus
    
    def read_arg = '--nano-raw'
    if (read_type == 'nano-corr') {
        read_arg = '--nano-corr'
    } else if (read_type == 'nano-hq') {
        read_arg = '--nano-hq'
    } else if (read_type == 'pacbio-corr') {
        read_type = '--pacbio-corr'
    } else {
        read_type = '--pacbio-corr'
    }
    
    """
    set -euxo pipefail
    
    flye ${read_arg} ${reads} --threads ${num_core} --out-dir asm
    
    mv asm/assembly.fasta ${prefix}.flye.fa
    mv asm/assembly_graph.gfa ${prefix}.flye.gfa
    """
}

params.read_type = "nano-raw"
    
workflow {
    Assemble(file(params.reads), params.prefix, params.read_type)
    
    emit:
    fa = Assemble.out[0]
    gfa = Assemble.out[1]
}
