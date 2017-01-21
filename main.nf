#!/usr/bin/env nextflow
//directory = '/projects/b1059/data/fastq/WI/dna/processed/**/'

/*
    Filtering configuration
*/
min_depth = 3
site_list=Channel.fromPath("CB4856.20160408.sitelist.tsv.gz")
site_list_index=Channel.fromPath("CB4856.20160408.sitelist.tsv.gz.tbi")
concordance_script=Channel.fromPath("concordance.R")

/*
    Set these parameters in nextflow.config
*/
tmpdir = config.tmpdir
reference = config.reference
cores = config.cores
analysis_dir = config.analysis_root + "/${params.type}-analysis"
bam_dir = config.bam_dir
call_variant_cpus = config.call_variant_cpus

// Define contigs here!
contig_list = ["I", "II", "III", "IV", "V", "X", "MtDNA"]
contigs = Channel.from(contig_list)

println "Processing RIL Data"
println "Using Reference: ${reference}" 

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
        cp ${strain_set_file} ${analysis_dir}/strain_set.json
    """
}

/*
    Fastq alignment
*/

process perform_alignment {

    echo true

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
        val "${fq_pair_id}" into fq_pair_id_concordance_val
        set file("${fq_pair_id}.bam"), file("${fq_pair_id}.bam.bai") into fq_idx_stats_bam
        set file("${fq_pair_id}.bam"), file("${fq_pair_id}.bam.bai") into fq_stat_bams
        set SM, file("${fq_pair_id}.bam"), file("${fq_pair_id}.bam.bai") into fq_pair_id_concordance

    
    """
        bwa mem -t ${cores} -R '${RG}' ${reference} ${fq1} ${fq2} | \\
        sambamba view --nthreads=${cores} --sam-input --format=bam --with-header /dev/stdin | \\
        sambamba sort --nthreads=${cores} --show-progress --tmpdir=${tmpdir} --out=${fq_pair_id}.bam /dev/stdin
        sambamba index --nthreads=${cores} ${fq_pair_id}.bam
    """
}

/* 
    Call variants at the individual level for concordance
*/

process fq_call_variants {

    tag { fq_pair_id }

    publishDir analysis_dir + "/fq_pair_vcf", mode: 'copy'

    input:
        val fq_pair_id from fq_pair_id_concordance_val
        set val(SM), file("${fq_pair_id}.bam"), file("${fq_pair_id}.bam.bai") from fq_pair_id_concordance

    output:
        file("out.tsv") into fq_individual_sites

    """
    # Perform individual-level calling
    contigs="`samtools view -H ${fq_pair_id}.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40`"
    echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${cores} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags DP,AD,ADF,ADR,SP --fasta-ref ${reference} ${fq_pair_id}.bam | bcftools call --skip-variants indels --variants-only --multiallelic-caller -O z  -  > ${fq_pair_id}.{}.individual.vcf.gz"
    order=`echo \${contigs} | tr ' ' '\\n' | awk '{ print "${fq_pair_id}." \$1 ".individual.vcf.gz" }'`
    
    # Output variant sites
    bcftools concat \${order} -O v | vk geno het-polarization - | bcftools view -O z > ${fq_pair_id}.individual.vcf.gz
    bcftools index ${fq_pair_id}.individual.vcf.gz
    rm \${order}

    bcftools view -M 2 -m 2 -O v ${fq_pair_id}.individual.vcf.gz | \\
    bcftools filter --include 'DP > 3' | \\
    grep -v "0\\/1" | \\
    bcftools query -f '%CHROM-%POS[\\t%GT]\\n' | \\
    awk '{ sub(/0\\/0/, "0");  sub(/1\\/1/, "1"); print "${SM}.${fq_pair_id}\\t" \$0 }' > out.tsv

    """
}

process fq_SM_concordance {

    publishDir analysis_dir + "/concordance", mode: 'copy'

    input:
        file("out?.tsv") from fq_individual_sites.toList()
        file(s:"script.R") from concordance_script

    output:
        file("fq_concordance.svg")
        file("fq_concordance.png")
        file("fq_concordance.tsv")

    """
        Rscript --vanilla ${s} 
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
        file("fq_bam_idxstats.tsv")

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > fq_bam_idxstats.tsv
        cat ${bam_idxstats.join(" ")} >> fq_bam_idxstats.tsv
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
        file("fq_bam_stats.tsv")

    """
        echo -e "fq_pair_id\\tvariable\\tvalue\\tcomment" > fq_bam_stats.tsv
        cat ${stat_files.join(" ")} >> fq_bam_stats.tsv
    """
}

/*
    Fastq coverage
*/
process fq_coverage {

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


process fq_coverage_merge {

    publishDir analysis_dir + "/fq", mode: 'copy'

    input:
        val fq_set from fq_coverage.toSortedList()

    output:
        file("fq_coverage.full.tsv")
        file("fq_coverage.tsv")

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > fq_coverage.full.tsv
        cat ${fq_set.join(" ")} >> fq_coverage.full.tsv

        cat <(echo -e 'fq\\tcoverage') <( cat fq_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > fq_coverage.tsv
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
        file 'SM_bam_idxstats' into bam_idxstats_set

    """
        samtools idxstats ${SM}.bam | awk '{ print "${SM}\\t" \$0 }' > SM_bam_idxstats
    """
}

process combine_idx_stats {

    publishDir analysis_dir +"/SM", mode: 'copy'

    input:
        val bam_idxstats from bam_idxstats_set.toSortedList()

    output:
        file("SM_bam_idxstats.tsv")

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > SM_bam_idxstats.tsv
        cat ${bam_idxstats.join(" ")} >> SM_bam_idxstats.tsv
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
        file("SM_bam_stats.tsv")

    """
        echo -e "fq_pair_id\\tvariable\\tvalue\\tcomment" > SM_bam_stats.tsv
        cat ${stat_files.join(" ")} >> SM_bam_stats.tsv
    """
}



process format_duplicates {

    publishDir analysis_dir + "/duplicates", mode: 'copy'

    input:
        val duplicates_set from duplicates_file.toSortedList()

    output:
        file("bam_duplicates.tsv")


    """
        echo -e 'filename\\tlibrary\\tunpaired_reads_examined\\tread_pairs_examined\\tsecondary_or_supplementary_rds\\tunmapped_reads\\tunpaired_read_duplicates\\tread_pair_duplicates\\tread_pair_optical_duplicates\\tpercent_duplication\\testimated_library_size' > bam_duplicates.tsv
        for i in ${duplicates_set.join(" ")}; do
            f=\$(basename \${i})
            cat \${i} | awk -v f=\${f/.duplicates.txt/} 'NR >= 8 && \$0 !~ "##.*" && \$0 != ""  { print f "\\t" \$0 } NR >= 8 && \$0 ~ "##.*" { exit }'  >> bam_duplicates.tsv
        done;
    """
}


/*
    Coverage Bam
*/
process SM_coverage {

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


process SM_coverage_merge {

    publishDir analysis_dir + "/SM", mode: 'copy'


    input:
        val sm_set from SM_coverage.toSortedList()

    output:
        file("SM_coverage.full.tsv")
        file("SM_coverage.tsv")

    """
        echo -e 'SM\\tcontig\\tstart\\tend\\tproperty\\tvalue' > SM_coverage.full.tsv
        cat ${sm_set.join(" ")} >> SM_coverage.full.tsv

        # Generate condensed version
        cat <(echo -e 'strain\\tcoverage') <(cat SM_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > SM_coverage.tsv
    """

}


/* 
    Call variants using the merged site list
*/

site_list =  site_list.spread(site_list_index)
union_vcf_channel = merged_bams_union.spread(site_list)



process call_variants_union {

    echo true

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
       file("union_vcfs.txt") into union_vcfs

    """
        echo ${vcf_set.join(" ")} | tr ' ' '\\n' > union_vcfs.txt
    """
}

union_vcfs_in = union_vcfs.spread(contigs)

process merge_union_vcf_chromosome {

    tag { chrom }

    input:
        set file(union_vcfs:"union_vcfs.txt"), val(chrom) from union_vcfs_in

    output:
        val(chrom) into contigs_list_in
        file("${chrom}.merged.raw.vcf.gz") into raw_vcf

    """
        bcftools merge --threads 10 --regions ${chrom} -O z -m all --file-list ${union_vcfs} > ${chrom}.merged.raw.vcf.gz
        bcftools index ${chrom}.merged.raw.vcf.gz
    """
}

// Generate a list of ordered files.
contig_raw_vcf = contig_list*.concat(".merged.raw.vcf.gz")

process concatenate_union_vcf {

    input:
        val merge_vcf from raw_vcf.toList()

    output:
        set file("merged.raw.vcf.gz"), file("merged.raw.vcf.gz.csi") into raw_vcf_concatenated

    """
        for i in ${merge_vcf.join(" ")}; do
            ln  -s \${i} `basename \${i}`;
        done;
        chrom_set="";
        bcftools concat -O z ${contig_raw_vcf.join(" ")}  > merged.raw.vcf.gz
        bcftools index merged.raw.vcf.gz
    """
}


process filter_union_vcf {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        set file("merged.raw.vcf.gz"), file("merged.raw.vcf.gz.csi") from raw_vcf_concatenated

    output:
        set file("merged.filtered.vcf.gz"), file("merged.filtered.vcf.gz.csi") into filtered_vcf

    """
        min_depth=${min_depth}

        bcftools view merged.raw.vcf.gz | \\
        vk geno het-polarization - | \\
        bcftools filter -O u --threads 16 --set-GTs . --include "FORMAT/DP > \${min_depth}" | \\
        bcftools view -O z - > merged.filtered.vcf.gz
        bcftools index -f merged.filtered.vcf.gz
    """
}

filtered_vcf.into { filtered_vcf_gtcheck; filtered_vcf_stat }

process gtcheck_tsv {

    publishDir analysis_dir + "/concordance", mode: 'copy'

    input:
        set file("merged.filtered.vcf.gz"), file("merged.filtered.vcf.gz.csi") from filtered_vcf_gtcheck

    output:
        file("SM.gtcheck.tsv") into gtcheck

    """
        echo -e "discordance\\tsites\\tavg_min_depth\\ti\\tj" > SM.gtcheck.tsv
        bcftools gtcheck -H -G 1 merged.filtered.vcf.gz | egrep '^CN' | cut -f 2-6 >> SM.gtcheck.tsv
    """

}


process stat_tsv {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        set file("merged.filtered.vcf.gz"), file("merged.filtered.vcf.gz.csi")  from filtered_vcf_stat

    output:
        file("filtered.stats.txt")

    """
        bcftools stats --verbose merged.filtered.vcf.gz > filtered.stats.txt
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
