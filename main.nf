#!/usr/bin/env nextflow
params.directory = '/projects/b1059/data/fastq/WI/dna/processed/**/'
params.analysis_dir = "/projects/b1059/analysis/WI_concordance"
params.kmer_size = 31
params.sketches = 20000

output_base = "k${params.kmer_size}_s${params.sketches}.tsv"
out_fq_tsv = "fq_${output_base}"
out_strain_tsv = "strain_${output_base}"
out_isotype_tsv = "isotype_${output_base}"

println "Running Concordance on " + params.directory
println "kmers: " + params.kmer_size
println "sketches: " + params.sketches

Channel.fromFilePairs(params.directory + '*{1,2}P.fq.gz', flat: true)
        .into { fq_pairs }

// Construct strain and isotype lists
import groovy.json.JsonSlurper
def strain_keys = []
def strain_values = []
def isotype_keys = []
def isotype_values = []

// Strain
def strainFile = new File('strain_set.json')
def strainJSON = new JsonSlurper().parseText(strainFile.text)
strainJSON.each{ key, value -> strain_keys << key; }
strainJSON.each{ key, value -> strain_values << value; }

// Isotype
def isoFile = new File('isotype_set.json')
def isoJSON = new JsonSlurper().parseText(isoFile.text)
isoJSON.each{ key, value -> isotype_keys << key; }
isoJSON.each{ key, value -> isotype_values << value; }


/*
    Fastq concordance
*/
process sketch_fq_files {

    cpus 16

    input:
        set dataset_id, file(fq1), file(fq2) from fq_pairs
    output:
        set file("${dataset_id}.msh") into sketches

    """
    zcat ${fq1} ${fq2} | pigz > ${dataset_id}.fq.gz
    mash sketch -r -p 16 -m 2 -k ${params.kmer_size} -s ${params.sketches} -o ${dataset_id} ${dataset_id}.fq.gz
    rm ${dataset_id}.fq.gz
    """
}

process combine_fq_sketch_files {

    input:
    file sketch from sketches.toList()

    output:
    file "output.msh" into output
    file "${out_tsv}"

    publishDir "/projects/b1059/analysis/WI_concordance/fq", mode: 'copy'

    """
    mash paste output ${sketch}
    mash dist output.msh output.msh > ${out_fq_tsv}
    """
}

/*
    Strain Concordance
*/
process sketch_strain_files {

    cpus 16

    input:
        val key from strain_keys
        file strain_set from strain_values
    output:
        set file("${key}.msh") into strain_sketches

    """
    zcat ${strain_set} | pigz > ${key}.fq.gz
    mash sketch -r -p 16 -m 2 -k ${params.kmer_size} -s ${params.sketches} -o ${key} ${key}.fq.gz
    rm ${key}.fq.gz
    """
}

process combine_strain_files {

    input:
    file sketch from strain_sketches.toList()

    output:
    file "output.msh" into output
    file "${out_strain_tsv}"

    publishDir "/projects/b1059/analysis/WI_concordance/strain", mode: 'copy'

    """
    mash paste output ${sketch}
    mash dist output.msh output.msh > ${out_strain_tsv}
    """
}


/*
    Isotype Concordance
*/
process sketch_isotype_files {

    cpus 16

    input:
        val key from isotype_keys
        file strain_set from isotype_values
    output:
        set file("${key}.msh") into isotype_sketches

    """
    zcat ${strain_set} | pigz > ${key}.fq.gz
    mash sketch -r -p 16 -m 2 -k ${params.kmer_size} -s ${params.sketches} -o ${key} ${key}.fq.gz
    rm ${key}.fq.gz
    """
}

process combine_isotype_files {

    input:
    file sketch from isotype_sketches.toList()

    output:
    file "output.msh" into output
    file "${out_strain_tsv}"

    publishDir "/projects/b1059/analysis/WI_concordance/isotype", mode: 'copy'

    """
    mash paste output ${sketch}
    mash dist output.msh output.msh > ${out_isotype_tsv}
    """
}


workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

