#!/usr/bin/env nextflow

params.directory = "$PWD/"
params.out = params.directory.replace("raw", "processed")
println params.out
params.threads = 8
println "Running Trimmomatic on " + params.directory
println params.directory + '*_R{1,2}_.*.fastq.gz'
Channel.fromFilePairs(params.directory + 'EA2_NIC277*_R{1,2}_001.fastq.gz', flat: true)
        .into { trimmomatic_read_pairs }


process make_out_dir {
    
    executor 'local'

    """
    mkdir -p ${params.out}
    """
}

process trim {

    publishDir params.out, mode: 'copy'

    cpus 8

    input:
        set dataset_id, file(forward), file(reverse) from trimmomatic_read_pairs

    output:
        set file("${dataset_id}_1P.fq.gz"), file("${dataset_id}_2P.fq.gz") into trim_output

    """
    trimmomatic PE -threads ${params.threads} $forward $reverse -baseout ${dataset_id}.fq.gz ILLUMINACLIP:/home/dec211/.linuxbrew/share/trimmomatic/adapters/NexteraPE-PE.fa:2:80:10 MINLEN:45
    rm ${dataset_id}_1U.fq.gz
    rm ${dataset_id}_2U.fq.gz
    """

}

process copy_fq_files {

    executor 'local'

    """
    cp ${params.directory}/.description ${params.out}/.description
    cp ${params.directory}/.fqdata ${params.out}/.fqdata
    """

}

process perform_fq_profile {
    input:
        set f1, f2 from trim_output

    """
    fq profile ${f1} ${f2}
    """
}

