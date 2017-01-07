#!/usr/bin/env nextflow
directory = '/projects/b1059/data/fastq/WI/dna/processed/**/'
tmpdir = "/projects/b1059/data/tmp"
genome="WS245"
reference = "/projects/b1059/data/genomes/c_elegans/${genome}/${genome}.fa.gz"
cores = 16
compression_threads = 4
date = new Date().format( 'yyyy-MM-dd' )
analysis_dir = "/projects/b1059/analysis/${date}-WI-summary"

println "Running Concordance on " + directory
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


process setup_dirs {

    executor 'local'

    """
        mkdir -p ${analysis_dir}
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

    publishDir analysis_dir, mode: 'move'

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
        val SM into merged_bam_SM
        val SM into merged_bam_SM_for_coverage
        set file("${SM}.bam"), file("${SM}.bam.bai") into merged_bams
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
        val SM from merged_bam_SM_for_coverage
        set file("${SM}.bam"), file("${SM}.bam.bai") from merged_bams_for_coverage

    output:
        file("${SM}.coverage.tsv") into SM_coverage


    """
        bam coverage ${SM}.bam > ${SM}.coverage.tsv
    """
}


process coverage_SM_merge {

    publishDir analysis_dir, mode: 'move'

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
    echo true
    input:
        val SM from merged_bam_SM
        set file("${SM}.bam"), file("${SM}.bam.bai") from merged_bams

    output:
        file("${SM}.individual.sites.tsv") into individual_sites

    """
    # Perform individual-level calling
    contigs="`samtools view -H ${SM}.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40`"
    echo \${contigs}
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

    publishDir analysisDir, mode: 'copy'
    
    input:
        val sites from individual_sites.toList()

    output:
        file("${date}.sitelist.tsv") 

    """
    cat ${sites.join(" ")} | sort -k1,1 -k2,2n | uniq > ${date}.sitelist.tsv
    """
}



workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
