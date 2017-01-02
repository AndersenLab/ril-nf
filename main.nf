#!/usr/bin/env nextflow

params.directory = '/projects/b1059/data/fastq/WI/dna/processed/**/'
params.analysis_dir = "/projects/b1059/analysis/WI_concordance"


println "Running Concordance on " + params.directory

Channel.fromFilePairs(params.directory + '*{1,2}P.fq.gz', flat: true)
        .into { fq_pairs }


process sketch_files {

    cpus 16

    input:
        set dataset_id, file(fq1), file(fq2) from fq_pairs
    output:
        set file("${dataset_id}.msh") into sketches

    """
    zcat ${fq1} ${fq2} > ${dataset_id}.fq.gz
    mash sketch -r -p 16 -m 2 -k 21 -s 20000 -o ${dataset_id} ${dataset_id}.fq.gz
    rm ${dataset_id}.fq.gz
    """
}

process combine_sketch_files {

    input:
    file sketch from sketches.toList()

    output:
    file "output.msh" into output
    file "concordance.tsv"

    publishDir "/projects/b1059/analysis/WI_concordance/sketches/fq", mode: 'copy'

    """
    mash paste output ${sketch}
    mash dist output.msh output.msh > concordance.tsv
    """
}


workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

