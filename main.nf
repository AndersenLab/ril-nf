#!/usr/bin/env nextflow
//directory = '/projects/b1059/data/fastq/WI/dna/processed/**/'

/*
    Filtering configuration
*/
min_depth = 3
site_list = Channel.fromPath("sitelist.tsv.gz")
site_list_index = Channel.fromPath("sitelist.tsv.gz.tbi")

/*
    Set these parameters in nextflow.config
/*
tmpdir = config.tmpdir
reference = config.reference
cores = config.cores
compression_threads = config.compression_threads
date = config.date
genome = config.genome
analysis_dir = config.analysis_dir

println "Processing RIL Data"
println "Using Reference: ${genome}" 

// Construct strain and isotype lists
import groovy.json.JsonSlurper

def strain_set = []

// Strain
def strainFile = new File('strain_set.json')
def strainJSON = new JsonSlurper().parseText(strainFile.text)

strainJSON.each { SM, RG ->
    RG.each { k, v ->
        strain_set << [SM, k, v[0], v[1], v[2]]
    }
}

strain_set_file = Channel.fromPath('strain_set.json')

process setup_dirs {

    executor 'local'

    input:
        file strain_set_file

    """
        mkdir -p ${analysis_dir}
        cp ${strain_set_file} ${analysis_dir}/${date}.strain_set.json
    """
}

/*
    Fastq alignment
*/

process perform_alignment {

    cpus 4

    tag { fq_pair_id }

    input:
        set SM, RG, fq1, fq2, fq_pair_id from strain_set
    output:
        set SM, file("${fq_pair_id}.bam"), file("${fq_pair_id}.bam.bai") into sample_aligned_bams
        val "${fq_pair_id}" into fq_pair_id_cov
        file "${fq_pair_id}.bam" into fq_cov_bam
        file "${fq_pair_id}.bam.bai" into fq_cov_bam_indices

    
    """
        bwa mem -t ${cores} -R '${RG}' ${reference} ${fq1} ${fq2} | \\
        sambamba view --nthreads=${cores} --sam-input --format=bam --with-header /dev/stdin | \\
        sambamba sort --nthreads=${cores} --show-progress --tmpdir=${tmpdir} --out=${fq_pair_id}.bam /dev/stdin
        sambamba index --nthreads=${cores} ${fq_pair_id}.bam
    """
}

/*
    Fastq coverage
*/
process coverage_fq {

    tag { fq_pair_id }

    input:
        val fq_pair_id from fq_pair_id_cov
        file("${fq_pair_id}.bam") from fq_cov_bam
        file("${fq_pair_id}.bam.bai") from fq_cov_bam_indices
    output:
        file("${fq_pair_id}.coverage.tsv") into fq_coverage


    """
        bam coverage ${fq_pair_id}.bam > ${fq_pair_id}.coverage.tsv
    """
}


process coverage_fq_merge {

    publishDir analysis_dir, mode: 'copy'

    input:
        val fq_set from fq_coverage.toSortedList()

    output:
        file("${date}.fq_coverage.full.tsv")
        file("${date}.fq_coverage.tsv")

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > ${date}.fq_coverage.tsv
        cat ${fq_set.join(" ")} >> ${date}.fq_coverage.full.tsv

        cat <(echo -e 'fq\\tcoverage') <( cat ${date}.fq_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > ${date}.fq_coverage.tsv
    """
}


process merge_bam {

    cpus cores

    tag { SM }

    input:
        set SM, bam, index from sample_aligned_bams.groupTuple()

    output:
        val SM into merged_SM_coverage
        val SM into merged_SM_individual
        val SM into merged_SM_union
        set file("${SM}.bam"), file("${SM}.bam.bai") into merged_bams_for_coverage
        set file("${SM}.bam"), file("${SM}.bam.bai") into merged_bams_union

    """

    count=`echo ${bam.join(" ")} | tr ' ' '\\n' | wc -l`

    if [ "\${count}" -eq "1" ]; then
        ln -s ${bam.join(" ")} ${SM}.merged.bam
        ln -s ${bam.join(" ")}.bai ${SM}.merged.bam.bai
    else
        sambamba merge --nthreads=${cores} --show-progress ${SM}.merged.bam ${bam.join(" ")}
        sambamba index --nthreads=${cores} ${SM}.merged.bam
    fi

    picard MarkDuplicates I=${SM}.merged.bam O=${SM}.bam M=${SM}.duplicates.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false
    sambamba index --nthreads=${cores} ${SM}.bam
    """
}


/*
    Coverage Bam
*/
process coverage_SM {

    tag { SM }

    input:
        val SM from merged_SM_coverage
        set file("${SM}.bam"), file("${SM}.bam.bai") from merged_bams_for_coverage

    output:
        val SM into SM_coverage_sample
        file("${SM}.coverage.tsv") into SM_coverage


    """
        bam coverage ${SM}.bam > ${SM}.coverage.tsv
    """
}


process coverage_SM_merge {

    publishDir analysis_dir, mode: 'copy'


    input:
        val sm_set from SM_coverage.toSortedList()

    output:
        file("${date}.SM_coverage.full.tsv")
        file("${date}.SM_coverage.tsv")

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > ${date}.SM_coverage.tsv
        cat ${sm_set.join(" ")} >> ${date}.SM_coverage.full.tsv

        # Generate condensed version
        cat <(echo -e 'strain\\tcoverage') <(cat ${date}.SM_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > ${date}.SM_coverage.tsv
    """

}


/* 
    Call variants using the merged site list
*/


process call_variants_union {

    cpus 6

    tag { SM }

    input:
        val SM from merged_SM_union
        set file("${SM}.bam"), file("${SM}.bam.bai") from merged_bams_union
        file 'sitelist.tsv.gz' from site_list 
        file 'sitelist.tsv.gz.tbi' from site_list_index

    output:
        val SM into union_vcf_SM
        file("${SM}.union.vcf.gz") into union_vcf_set
        file("${SM}.union.vcf.gz.csi") into union_vcf_set_indices


    """
        contigs="`samtools view -H ${SM}.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40`"
        echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${cores} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags DP,AD,ADF,ADR,INFO/AD,SP --fasta-ref ${reference} ${SM}.bam | bcftools call -T sitelist.tsv.gz --skip-variants indels --variants-only --multiallelic-caller -O z  -  > ${SM}.{}.union.vcf.gz"
        order=`echo \${contigs} | tr ' ' '\\n' | awk '{ print "${SM}." \$1 ".union.vcf.gz" }'`

        # Output variant sites
        bcftools concat \${order} -O v | vk geno het-polarization - | bcftools view -O z > ${SM}.union.vcf.gz
        bcftools index ${SM}.union.vcf.gz
        rm \${order}
    """

}


process generate_union_vcf_list {

    cpus 1 

    publishDir analysis_dir, mode: 'copy'

    input:
       val vcf_set from union_vcf_set.toSortedList()

    output:
       file("${date}.union_vcfs.txt") into union_vcfs

    """
        echo ${vcf_set.join(" ")} | tr ' ' '\\n' > ${date}.union_vcfs.txt
    """
}


process merge_union_vcf {

    cpus cores

    publishDir analysis_dir, mode: 'copy'

    input:
        val SM from union_vcf_SM.toSortedList()
        file(union_vcfs:"union_vcfs.txt") from union_vcfs

    output:
        file("${date}.merged.raw.vcf.gz") into raw_vcf
        file("${date}.merged.filtered.vcf.gz") into filtered_vcf

    """
        bcftools merge --threads 24 -O z -m all --file-list ${union_vcfs} > ${date}.merged.raw.vcf.gz
        bcftools index ${date}.merged.raw.vcf.gz

        min_depth=${min_depth}
        bcftools view ${date}.merged.raw.vcf.gz | \\
        bcftools filter -O u --threads 16 --set-GTs . --include "FORMAT/DP > \${min_depth}" | \\
        bcftools view -O z - > ${date}.merged.filtered.vcf.gz
        bcftools index -f ${date}.merged.filtered.vcf.gz

    """

}

filtered_vcf.into { filtered_vcf_gtcheck; filtered_vcf_stat }

process gtcheck_tsv {

    publishDir analysis_dir, mode: 'copy'

    input:
        file("${date}.merged.filtered.vcf.gz") from filtered_vcf_gtcheck

    output:
        file("${date}.gtcheck.tsv") into gtcheck

    """
        echo -e "discordance\\tsites\\tavg_min_depth\\ti\\tj" > ${date}.gtcheck.tsv
        bcftools gtcheck -H -G 1 ${date}.merged.filtered.vcf.gz | egrep '^CN' | cut -f 2-6 >> ${date}.gtcheck.tsv
    """

}


process stat_tsv {

    publishDir analysis_dir, mode: 'copy'

    input:
        file("${date}.merged.filtered.vcf.gz") from filtered_vcf_stat

    output:
        file("${date}.stats.txt")

    """
        bcftools stats ${date}.merged.filtered.vcf.gz > ${date}.stats.txt
    """

}


workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
