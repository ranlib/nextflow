include { Minimap2 } from './align_reads.nf'

params.tags_to_preserve = []
params.library = ''

workflow {
    Minimap2(params.reads, params.ref_fasta, params.RG, params.map_preset, params.library, params.tags_to_preserve, params.prefix)
}
