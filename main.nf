#!/usr/bin/env nextflow

params.directory = "$PWD/"
println "Running Fastq Profiler on " + params.directory

fastqs = Channel.fromPath( params.directory + '*.fastq.gz' )

process profile_fastqs {

    input:
    file 'query' from fastqs

    """
    fq profile --fastqc ${query} 
    """
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
