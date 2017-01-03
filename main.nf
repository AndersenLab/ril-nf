#!/usr/bin/env nextflow
directory = '/projects/b1059/data/fastq/WI/dna/processed/**/'
analysis_dir = "/projects/b1059/analysis/WI_concordance"
kmer_size = 21
sketch_count = 20000
min_kmer_count = 4

output_base = "k${kmer_size}_s${sketch_count}_m${min_kmer_count}.tsv"
out_fq_tsv = "fq_${output_base}"
out_strain_tsv = "strain_${output_base}"
out_isotype_tsv = "isotype_${output_base}"

println "Running Concordance on " + directory
println "kmers: " + kmer_size
println "sketches: " + sketch_count
println "minimum count of kmer: " + min_kmer_count

Channel.fromFilePairs(directory + '*{1,2}P.fq.gz', flat: true)
        .into { fq_pairs }

// Construct strain and isotype lists
import groovy.json.JsonSlurper
def strain_set = []
def isotype_set = []

// Strain
def strainFile = new File('strain_set.json')
def strainJSON = new JsonSlurper().parseText(strainFile.text)

strainJSON.each { k, v -> for (i in v) {
    strain_set << [k, i]
    }
}

// Isotype
def isoFile = new File('isotype_set.json')
def isoJSON = new JsonSlurper().parseText(isoFile.text)

isoJSON.each { k, v -> for (i in v) {
    isotype_set << [k, i]
    }
}

strain_ch = Channel.from(strain_set)
       .groupTuple()

isotype_ch = Channel.from(isotype_set)
       .groupTuple()


/*
    Fastq concordance
*/
process sketch_fq_files {

    cpus 16

    input:
        set dataset_id, file(fq1), file(fq2) from fq_pairs
    output:
        file("${dataset_id}.msh") into fq_sketches

    """
    zcat ${fq1} ${fq2} | pigz > ${dataset_id}.fq.gz
    mash sketch -r -p 16 - ${min_kmer_count} -k ${kmer_size} -s ${sketch_count} -o ${dataset_id} ${dataset_id}.fq.gz
    rm ${dataset_id}.fq.gz
    """
}

process combine_fq_sketch_files {

    input:
    file sketch from fq_sketches.toList()

    output:
    file "output.msh" into fq_output
    file "${out_fq_tsv}"

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
    echo true
    cpus 16

    input:
        set seq_id, fq from strain_ch
    output:
        file("${seq_id}.msh") into strain_sketches

    """
    zcat ${fq.join(" ")} | pigz > ${seq_id}.fq.gz
    mash sketch -r -p 16 - ${min_kmer_count} -k ${kmer_size} -s ${sketch_count} -o ${seq_id} ${seq_id}.fq.gz
    rm ${seq_id}.fq.gz
    """
}

process combine_strain_files {

    input:
    file sketch from strain_sketches.toList()

    output:
    file "output.msh" into sketch_output
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
        set seq_id, fq from isotype_ch
    output:
        file("${seq_id}.msh") into isotype_sketches

    """
    zcat ${fq.join(" ")} | pigz > ${seq_id}.fq.gz
    mash sketch -r -p 16 - ${min_kmer_count} -k ${kmer_size} -s ${sketch_count} -o ${seq_id} ${seq_id}.fq.gz
    rm ${seq_id}.fq.gz
    """
}

process combine_isotype_files {

    input:
    file sketch from isotype_sketches.toList()

    output:
    file "output.msh" into isotype_output
    file "${out_isotype_tsv}"

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

