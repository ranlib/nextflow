params.SRRs = "./*.fq.gz"

workflow {

  Channel
    .fromPath( params.SRRs )
    .view()
}