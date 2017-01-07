#!/usr/bin/env nextflow
//directory = '/projects/b1059/data/fastq/WI/dna/processed/**/'

/*
    Filtering configuration
*/
min_depth=3
qual=30
mq=40
dv_dp=0.5

tmpdir = "/projects/b1059/data/tmp"
genome = "WS245"
reference = "/projects/b1059/data/genomes/c_elegans/${genome}/${genome}.fa.gz"
cores = 16
compression_threads = 4
date = new Date().format( 'yyyy-MM-dd' )
analysis_dir = "/projects/b1059/analysis/${date}-WI-summary"

tmpdir = config.tmpdir
reference = config.reference
cores = config.cores
compression_threads = config.compression_threads
date = config.date
genome = config.genome
analysis_dir = config.analysis_dir

println "Running Concordance on Wild Isolates"
println "Using Reference: ${genome}" 

// Construct strain and isotype lists
import groovy.json.JsonSlurper

def strain_set = []
def isotype_set = []

// Strain
def strainFile = new File('strain_set.json')
def strainJSON = new JsonSlurper().parseText(strainFile.text)

strainJSON.each { SM, RG ->
    RG.each { k, v ->
        strain_set << [SM, k, v[0], v[1], v[2]]
    }
}

strain_set = strain_set[1..15]

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
    Fastq concordance
*/
process perform_alignment {

    cpus cores

    // load gcc module so pthread_cancel works!
    module 'gcc/5.1.0'

    input:
        set SM, RG, fq1, fq2, fq_pair_id from strain_set
    output:
        set SM, file("${fq_pair_id}.bam"), file("${fq_pair_id}.bam.bai") into sample_aligned_bams
        val "${fq_pair_id}" into fq_pair_id_cov
        file "${fq_pair_id}.bam" into fq_cov_bam
        file "${fq_pair_id}.bam.bai" into fq_cov_bam_indices

    
    """
        bwa mem -t ${cores} -R '${RG}' ${reference} ${fq1} ${fq2} | \
        head -n 200000 | \
        sambamba view --nthreads=${cores} --sam-input --format=bam --with-header /dev/stdin | \
        sambamba sort --nthreads=${cores} --show-progress --tmpdir=${tmpdir} --out=${fq_pair_id}.bam /dev/stdin
        sambamba index --nthreads=${cores} ${fq_pair_id}.bam
    """
}

/*
    Fastq coverage
*/
process coverage_fq {

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
        val fq_set from fq_coverage.toList()

    output:
        file("${date}.fq_coverage.tsv")

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > ${date}.fq_coverage.tsv
        cat ${fq_set.join(" ")} >> ${date}.fq_coverage.tsv
    """
}


process merge_bam {
    echo true

    cpus cores


    // load gcc module so pthread_cancel works!
    module 'gcc/5.1.0'

    input:
        set SM, bam, index from sample_aligned_bams.groupTuple()

    output:
        val SM into merged_SM_coverage
        val SM into merged_SM_individual
        val SM into merged_SM_union
        set file("${SM}.bam"), file("${SM}.bam.bai") into merged_bams_individual
        set file("${SM}.bam"), file("${SM}.bam.bai") into merged_bams_union
        set file("${SM}.bam"), file("${SM}.bam.bai") into merged_bams_for_coverage

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

    input:
        val SM from merged_SM_coverage
        set file("${SM}.bam"), file("${SM}.bam.bai") from merged_bams_for_coverage

    output:
        file("${SM}.coverage.tsv") into SM_coverage


    """
        bam coverage ${SM}.bam > ${SM}.coverage.tsv
    """
}


process coverage_SM_merge {

    publishDir analysis_dir, mode: 'copy'

    input:
        val sm_set from SM_coverage.toList()

    output:
        file("${date}.SM_coverage.tsv")

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > ${date}.SM_coverage.tsv
        cat ${sm_set.join(" ")} >> ${date}.SM_coverage.tsv
    """

}


process call_variants_individual {

    cpus cores

    input:
        val SM from merged_SM_individual
        set file("${SM}.bam"), file("${SM}.bam.bai") from merged_bams_individual

    output:
        file("${SM}.individual.sites.tsv") into individual_sites

    """
    # Perform individual-level calling
    contigs="`samtools view -H ${SM}.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40`"
    echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${cores} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags AD,ADF,ADR,INFO/AD,SP --fasta-ref ${reference} ${SM}.bam | bcftools call --skip-variants indels --variants-only --multiallelic-caller -O z  -  > ${SM}.{}.individual.vcf.gz"
    order=`echo \${contigs} | tr ' ' '\\n' | awk '{ print "${SM}." \$1 ".individual.vcf.gz" }'`
    
    # Output variant sites
    bcftools concat \${order} -O v | vk geno het-polarization - | bcftools view -O z > ${SM}.individual.vcf.gz
    bcftools index ${SM}.individual.vcf.gz
    rm \${order}

    bcftools view -M 2 -m 2 -O v ${SM}.individual.vcf.gz | \\
    bcftools filter --include 'DP > 3' | \\
    egrep '(^#|1/1)' | \\
    bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\n' > ${SM}.individual.sites.tsv

    """
}

/*
    Merge individual sites
*/

process merge_variant_list {

    publishDir analysis_dir, mode: 'copy'
    
    input:
        val sites from individual_sites.toList()

    output:
        file("${date}.sitelist.tsv") into sitelist
        set file("${date}.sitelist.tsv.gz"), file("${date}.sitelist.tsv.gz.tbi") into gz_sitelist

    """
        cat ${sites.join(" ")} | sort -k1,1 -k2,2n | uniq > ${date}.sitelist.tsv
        bgzip ${date}.sitelist.tsv -c > ${date}.sitelist.tsv.gz && tabix -s1 -b2 -e2 ${date}.sitelist.tsv.gz
    """
}

/* 
    Call variants using the merged site list
*/

process call_variants_union {

    cpus cores

    input:
        val SM from merged_SM_union
        set file("${SM}.bam"), file("${SM}.bam.bai") from merged_bams_union
        set file('sitelist.tsv.gz'), file('sitelist.tsv.gz.tbi') from gz_sitelist

    output:
        val SM into union_vcf_set_SM
        file("${SM}.union.vcf.gz") into union_vcf_set
        file("${SM}.union.vcf.gz.csi") into union_vcf_set_indices

    """
        contigs="`samtools view -H ${SM}.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40`"
        echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${cores} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags AD,ADF,ADR,INFO/AD,SP --fasta-ref ${reference} ${SM}.bam | bcftools call -T sitelist.tsv.gz --skip-variants indels --variants-only --multiallelic-caller -O z  -  > ${SM}.{}.union.vcf.gz"
        order=`echo \${contigs} | tr ' ' '\\n' | awk '{ print "${SM}." \$1 ".union.vcf.gz" }'`

        # Output variant sites
        bcftools concat \${order} -O v | vk geno het-polarization - | bcftools view -O z > ${SM}.union.vcf.gz
        bcftools index ${SM}.union.vcf.gz
        rm \${order}
    """

}


union_vcf_channel = Channel.create()

union_vcf_set.each { f ->
  union_vcf_channel.bind([f])
}


process merge_union_vcf {

    cpus cores

    publishDir analysis_dir, mode: 'copy'

    input:
        val SM from union_vcf_set_SM
        val union_vcfs from union_vcf_channel
        file("${SM}.union.vcf.gz.csi") from union_vcf_set_indices

    output:
        file("${date}.merged.raw.vcf.gz") into raw_vcf
        file("${date}.merged.filtered.vcf.gz") into filtered_vcf

    """
        bcftools merge --threads 24 -O z -m all ${union_vcfs.join(" ")} > ${date}.merged.no_filter.vcf.gz
        bcftools index ${date}.merged.raw.vcf.gz

        min_depth=${min_depth}
        qual=${qual}
        mq=${mq}
        dv_dp=${dv_dp}

        bcftools view ${date}.merged.raw.vcf.gz | \\
        bcftools filter -O u --threads 16 --set-GTs . --include "QUAL >= \${qual} || FORMAT/GT == '0/0'" |  \\
        bcftools filter -O u --threads 16 --set-GTs . --include "FORMAT/DP > \${min_depth}" | \\
        bcftools filter -O u --threads 16 --set-GTs . --include "INFO/MQ > \${mq}" | \\
        bcftools filter -O u --threads 16 --set-GTs . --include "(FORMAT/DV)/(FORMAT/DP) >= \${dv_dp} || FORMAT/GT == '0/0'" | \\
        vk filter REF --min=1 - | \\
        vk filter ALT --min=1 - | \\
        vk filter ALT --max=0.99 - | \\
        vk filter MISSING --max=0.90 --soft-filter='high_missing' --mode=x - | \\
        vk filter HET --max=0.10 --soft-filter='high_heterozygosity' --mode=+ - | \\
        bcftools view -O z - > ${date}.merged.filtered.vcf.gz
        bcftools index -f ${date}.merged.filtered.vcf.gz

    """

}

/*
process stat_vcf {

    input:
        file("${date}.merged.filtered.vcf.gz") from filtered_vcf

    output:

    """
    bcftools 

}
*/

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
