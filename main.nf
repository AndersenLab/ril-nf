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
*/
tmpdir = config.tmpdir
reference = config.reference
cores = config.cores
compression_threads = config.compression_threads
date = config.date
genome = config.genome
analysis_dir = config.analysis_dir
bam_dir = config.bam_dir
call_variant_cpus = config.call_variant_cpus

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

process workflow_info {

    executor 'local'

    publishDir analysis_dir, mode: 'copy'

    output:
        file("workflow.txt")

    """
        echo "Git info: $workflow.repository" > workflow.txt
        echo "Revision: $workflow.revision [$workflow.commitId]" >> workflow.txt
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
        val "${fq_pair_id}" into fq_stat_bam_val
        val "${fq_pair_id}" into fq_pair_id_idxstats
        set file("${fq_pair_id}.bam"), file("${fq_pair_id}.bam.bai") into fq_idx_stats_bam
        set file("${fq_pair_id}.bam"), file("${fq_pair_id}.bam.bai") into fq_stat_bams

    
    """
        bwa mem -t ${cores} -R '${RG}' ${reference} ${fq1} ${fq2} | \\
        sambamba view --nthreads=${cores} --sam-input --format=bam --with-header /dev/stdin | \\
        sambamba sort --nthreads=${cores} --show-progress --tmpdir=${tmpdir} --out=${fq_pair_id}.bam /dev/stdin
        sambamba index --nthreads=${cores} ${fq_pair_id}.bam
    """
}


/*
    fq idx stats
*/

process fq_idx_stats {
    
    input:
        val fq_pair_id from fq_pair_id_idxstats
        set file("${fq_pair_id}.bam"), file("${fq_pair_id}.bam.bai") from fq_idx_stats_bam
    output:
        file fq_idxstats into fq_idxstats_set

    """
        samtools idxstats ${fq_pair_id}.bam | awk '{ print "${fq_pair_id}\\t" \$0 }' > fq_idxstats
    """
}

process fq_combine_idx_stats {

    publishDir analysis_dir + "/fq", mode: 'copy'

    input:
        val bam_idxstats from fq_idxstats_set.toSortedList()

    output:
        file("${date}.fq_bam_idxstats.tsv")

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > ${date}.fq_bam_idxstats.tsv
        cat ${bam_idxstats.join(" ")} >> ${date}.fq_bam_idxstats.tsv
    """

}

/*
    fq bam stats
*/

process fq_bam_stats {

    tag { fq_pair_id }

    input:
        val fq_pair_id from fq_stat_bam_val
        set file("${fq_pair_id}.bam"), file("${fq_pair_id}.bam.bai") from fq_stat_bams

    output:
        file 'bam_stat' into bam_stat_files

    """
        cat <(samtools stats ${fq_pair_id}.bam | grep ^SN | cut -f 2- | awk '{ print "${fq_pair_id}\t" \$0 }' | sed 's/://g') > bam_stat
    """
}

process combine_bam_stats {

    publishDir analysis_dir + "/fq", mode: 'copy'

    input:
        val stat_files from bam_stat_files.toSortedList()

    output:
        file("${date}.fq_bam_stats.tsv")

    """
        echo -e "fq_pair_id\\tvariable\\tvalue\\tcomment" > ${date}.fq_bam_stats.tsv
        cat ${stat_files.join(" ")} >> ${date}.fq_bam_stats.tsv
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

    publishDir analysis_dir + "/fq", mode: 'copy'

    input:
        val fq_set from fq_coverage.toSortedList()

    output:
        file("${date}.fq_coverage.full.tsv")
        file("${date}.fq_coverage.tsv")

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > ${date}.fq_coverage.full.tsv
        cat ${fq_set.join(" ")} >> ${date}.fq_coverage.full.tsv

        cat <(echo -e 'fq\\tcoverage') <( cat ${date}.fq_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > ${date}.fq_coverage.tsv
    """
}


process merge_bam {

    publishDir bam_dir, mode: 'copy', pattern: '*.bam*'

    cpus cores

    tag { SM }

    input:
        set SM, bam, index from sample_aligned_bams.groupTuple()

    output:
        val SM into merged_SM_coverage
        val SM into merged_SM_individual
        val SM into merged_SM_union
        val SM into merged_SM_idxstats
        val SM into merged_SM_bamstats
        set file("${SM}.bam"), file("${SM}.bam.bai") into merged_bams_for_coverage
        set file("${SM}.bam"), file("${SM}.bam.bai") into merged_bams_union
        set file("${SM}.bam"), file("${SM}.bam.bai") into bams_idxstats
        set file("${SM}.bam"), file("${SM}.bam.bai") into bams_stats
        file("${SM}.duplicates.txt") into duplicates_file

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
 SM idx stats
*/

process idx_stats_SM {
    
    input:
        val SM from merged_SM_idxstats
        set file("${SM}.bam"), file("${SM}.bam.bai") from bams_idxstats
    output:
        file bam_idxstats into bam_idxstats_set

    """
        samtools idxstats ${SM}.bam | awk '{ print "${SM}\\t" \$0 }' > bam_idxstats
    """
}

process combine_idx_stats {

    publishDir analysis_dir +"/SM", mode: 'copy'

    input:
        val bam_idxstats from bam_idxstats_set.toSortedList()

    output:
        file("${date}.bam_idxstats.tsv")

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > ${date}.bam_idxstats.tsv
        cat ${bam_idxstats.join(" ")} >> ${date}.bam_idxstats.tsv
    """

}


/*
    SM bam stats
*/

process SM_bam_stats {

    tag { SM }

    input:
        val SM from merged_SM_bamstats
        set file("${SM}.bam"), file("${SM}.bam.bai") from bams_stats

    output:
        file 'bam_stat' into SM_bam_stat_files

    """
        cat <(samtools stats ${SM}.bam | grep ^SN | cut -f 2- | awk '{ print "${SM}\t" \$0 }' | sed 's/://g') > bam_stat
    """
}

process combine_SM_bam_stats {

    publishDir analysis_dir + "/SM", mode: 'copy'

    input:
        val stat_files from SM_bam_stat_files.toSortedList()

    output:
        file("${date}.SM_bam_stats.tsv")

    """
        echo -e "fq_pair_id\\tvariable\\tvalue\\tcomment" > ${date}.SM_bam_stats.tsv
        cat ${stat_files.join(" ")} >> ${date}.SM_bam_stats.tsv
    """
}



process format_duplicates {

    publishDir analysis_dir + "/duplicates", mode: 'copy'

    input:
        val duplicates_set from duplicates_file.toSortedList()

    output:
        file("${date}.bam_duplicates.tsv")


    """
        echo -e 'filename\\tlibrary\\tunpaired_reads_examined\\tread_pairs_examined\\tsecondary_or_supplementary_rds\\tunmapped_reads\\tunpaired_read_duplicates\\tread_pair_duplicates\\tread_pair_optical_duplicates\\tpercent_duplication\\testimated_library_size' > ${date}.bam_duplicates.tsv
        for i in ${duplicates_set.join(" ")}; do
            f=\$(basename \${i})
            cat \${i} | awk -v f=\${f/.duplicates.txt/} 'NR >= 8 && \$0 !~ "##.*" && \$0 != ""  { print f "\\t" \$0 } NR >= 8 && \$0 ~ "##.*" { exit }'  >> ${date}.bam_duplicates.tsv
        done;
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

    publishDir analysis_dir + "/SM", mode: 'copy'


    input:
        val sm_set from SM_coverage.toSortedList()

    output:
        file("${date}.SM_coverage.full.tsv")
        file("${date}.SM_coverage.tsv")

    """
        echo -e 'SM\\tcontig\\tstart\\tend\\tproperty\\tvalue' > ${date}.SM_coverage.full.tsv
        cat ${sm_set.join(" ")} >> ${date}.SM_coverage.full.tsv

        # Generate condensed version
        cat <(echo -e 'strain\\tcoverage') <(cat ${date}.SM_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > ${date}.SM_coverage.tsv
    """

}


/* 
    Call variants using the merged site list
*/

site_list =  site_list.spread(site_list_index)
union_vcf_channel = merged_bams_union.spread(site_list)



process call_variants_union {

    cpus call_variant_cpus

    tag { SM }

    input:
        val SM from merged_SM_union
        set file("${SM}.bam"), file("${SM}.bam.bai"), file('sitelist.tsv.gz'), file('sitelist.tsv.gz.tbi') from union_vcf_channel

    output:
        val SM into union_vcf_SM
        file("${SM}.union.vcf.gz") into union_vcf_set
        file("${SM}.union.vcf.gz.csi") into union_vcf_set_indices


    """
        contigs="`samtools view -H ${SM}.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40`"
        echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${cores} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags DP,AD,ADF,ADR,INFO/AD,SP --fasta-ref ${reference} ${SM}.bam | bcftools call -T sitelist.tsv.gz --skip-variants indels  --multiallelic-caller -O z  -  > ${SM}.{}.union.vcf.gz"
        order=`echo \${contigs} | tr ' ' '\\n' | awk '{ print "${SM}." \$1 ".union.vcf.gz" }'`

        # Output variant sites
        bcftools concat \${order} -O v | vk geno het-polarization - | bcftools view -O z > ${SM}.union.vcf.gz
        bcftools index ${SM}.union.vcf.gz
        rm \${order}
    """

}


process generate_union_vcf_list {

    cpus 1 

    publishDir analysis_dir + "/vcf", mode: 'copy'

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

    publishDir analysis_dir + "/vcf", mode: 'copy'

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

    publishDir analysis_dir + "/concordance", mode: 'copy'

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

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        file("${date}.merged.filtered.vcf.gz") from filtered_vcf_stat

    output:
        file("${date}.filtered.stats.txt")

    """
        bcftools stats --verbose ${date}.merged.filtered.vcf.gz > ${date}.filtered.stats.txt
    """

}


workflow.onComplete {
    def subject = 'RIL Workflow'
    def recipient = config.email

    ['mail', '-s', subject, recipient].execute() << """

    RIL Pipeline complete
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    profile: ${workflow.profile}
    Analysis Directory: ${analysis_dir}
    """
}
