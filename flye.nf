process Assemble {
    
    publishDir "flye", mode: 'copy'
    container "dbest/flye:v2.9.3"

    input:
    path reads
    val prefix
    val read_type
    val num_core
    
    output:
    path "${prefix}.flye.fa"
    path "${prefix}.flye.gfa"
    
    script:
    """
    set -euxo pipefail
    
    flye ${read_type} ${reads} --threads ${num_core} --out-dir asm
    
    mv asm/assembly.fasta ${prefix}.flye.fa
    mv asm/assembly_graph.gfa ${prefix}.flye.gfa
    """
}

params.read_type = "nano-raw"
    
workflow {
    Assemble(file(params.reads), params.prefix, params.read_type, params.num_core)
    fa = Assemble.out[0]
    gfa = Assemble.out[1]
}
