#!/usr/bin/env nextflow
//directory = '/projects/b1059/data/fastq/WI/dna/processed/**/'
/*
    Filtering configuration
*/
site_list=Channel.fromPath("CB4856.20160408.sitelist.tsv.gz")
site_list_index=Channel.fromPath("CB4856.20160408.sitelist.tsv.gz.tbi")
hmm_plot_script=Channel.fromPath("plot_hmm.R")
cross_object_script=file("generate_cross_object.R")

/*
    Set these parameters in nextflow.config
*/
date = new Date().format( 'yyyy-MM-dd' )
tmpdir = config.tmpdir
reference = config.reference
cores = config.cores
analysis_dir = config.analysis_dir
bam_dir = config.bam_dir
call_variant_cpus = config.call_variant_cpus

// Define contigs here!
contig_list = ["I", "II", "III", "IV", "V", "X", "MtDNA"]
contigs = Channel.from(contig_list)

println "Processing RIL Data"
println "Using Reference: ${reference}" 

// Construct strain and isotype lists
import groovy.json.JsonSlurper

strainFile = new File("fq_ril_sheet.tsv")
fqs = Channel.from(strainFile.collect { it.tokenize( '\t' ) })

process setup_dirs {

    executor 'local'

    publishDir analysis_dir, mode: 'copy'

    input:
        file("fq_ril_sheet.tsv") from Channel.fromPath("fq_ril_sheet.tsv")

    output:
        file("fq_ril_sheet.tsv")

    """
    """
}

/*
    Fastq alignment
*/

process perform_alignment {

    echo true

    cpus 4

    tag { ID }

    input:
        set SM, ID, LB, fq1, fq2 from fqs
    output:
        set SM, ID, LB, file("${ID}.bam"), file("${ID}.bam.bai") into aligned_bams
        set SM, file("${ID}.bam") into sample_aligned_bams
    
    """
        bwa mem -t ${cores} -R '@RG\tID:${ID}\tLB:${LB}\tSM:${SM}' ${reference} ${fq1} ${fq2} | \\
        sambamba view --nthreads=${cores} --sam-input --format=bam --with-header /dev/stdin | \\
        sambamba sort --nthreads=${cores} --show-progress --tmpdir=${tmpdir} --out=${ID}.bam /dev/stdin
        sambamba index --nthreads=${cores} ${ID}.bam
    """
}

aligned_bams.into { 
                           sample_bams_fq_idx_stats;
                           fq_stat_bams;
                           fq_cov_bam_indices; 
                         }

/*
    fq idx stats
*/

process fq_idx_stats {
    
    tag { ID }

    input:
        set SM, ID, LB, file("${ID}.bam"), file("${ID}.bam.bai") from sample_bams_fq_idx_stats
    output:
        file 'fq_idxstats' into fq_idxstats_set

    """
        samtools idxstats ${ID}.bam | awk '{ print "${ID}\\t" \$0 }' > fq_idxstats
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

    tag { ID }

    input:
        set SM, ID, LB, file("${ID}.bam"), file("${ID}.bam.bai") from fq_stat_bams

    output:
        file 'bam_stat' into bam_stat_files

    """
        cat <(samtools stats ${ID}.bam | grep ^SN | cut -f 2- | awk '{ print "${ID}\t" \$0 }' | sed 's/://g') > bam_stat
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

    tag { ID }

    input:
        set SM, ID, LB, file("${ID}.bam"), file("${ID}.bam.bai") from fq_cov_bam_indices
    output:
        file("${ID}.coverage.tsv") into fq_coverage


    """
        bam coverage ${ID}.bam > ${ID}.coverage.tsv
    """
}

process fq_coverage_merge {

    publishDir analysis_dir + "/fq", mode: 'copy'

    input:
        val fq_set from fq_coverage.toSortedList()

    output:
        file("fq_coverage.full.tsv")
        file("fq_coverage.tsv") into fq_coverage_plot

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > fq_coverage.full.tsv
        cat ${fq_set.join(" ")} >> fq_coverage.full.tsv

        cat <(echo -e 'fq\\tcoverage') <( cat fq_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > fq_coverage.tsv
    """
}

sample_aligned_bams.into { sample_aligned_bams_out; sample_aligned_bams_use}

sample_aligned_bams_out.groupTuple().println()

process merge_bam {

    echo true

    storeDir bam_dir

    cpus cores

    tag { SM }

    input:
        set SM, bam from sample_aligned_bams_use.groupTuple()

    output:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") into merged_SM
        file("${SM}.duplicates.txt") into duplicates_file

    """
    count=`echo ${bam.join(" ")} | tr ' ' '\\n' | wc -l`

    if [ "\${count}" -eq "1" ]; then
        ln -s ${bam[0]} ${SM}.merged.bam
        ln -s ${bam[0]}.bai ${SM}.merged.bam.bai
    else
        sambamba merge --nthreads=${cores} --show-progress ${SM}.merged.bam ${bam.join(" ")}
        sambamba index --nthreads=${cores} ${SM}.merged.bam
    fi

    picard MarkDuplicates I=${SM}.merged.bam O=${SM}.bam M=${SM}.duplicates.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false
    sambamba index --nthreads=${cores} ${SM}.bam
    """
}

merged_SM.into { 
                bams_stats;
                bams_idxstats;
                merged_bams_for_coverage;
                merged_bams_union;
                fq_concordance_bams;
                }



/*
 SM idx stats
*/

process idx_stats_SM {
    
    input:
        set SM, file("${SM}.bam"), file("${SM}.bam.bai") from bams_idxstats
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
        file("SM_bam_idxstats.tsv") into SM_bam_idxstats_plot

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
        set SM, file("${SM}.bam"), file("${SM}.bam.bai") from bams_stats

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
        file("bam_duplicates.tsv") into bam_duplicates_plot

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
        set SM, file("${SM}.bam"), file("${SM}.bam.bai") from merged_bams_for_coverage

    output:
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
        file("SM_coverage.tsv") into SM_coverage_plot
        file("SM_coverage.tsv") into SM_coverage_cross_obj

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
site_list.into { site_list_one; site_list_two }
union_vcf_channel = merged_bams_union.spread(site_list_one)
fq_concordance_sitelist = fq_concordance_bams.spread(site_list_two)


/* 
    Call variants at the individual level for concordance
*/

fq_concordance_script = file("fq_concordance.R")

process fq_concordance {

    cpus call_variant_cpus

    tag { SM }

    input:
        set val(SM), file("input.bam"), file("input.bam.bai"), file('sitelist.tsv.gz'), file('sitelist.tsv.gz.tbi') from fq_concordance_sitelist

    output:
        file('out.tsv') into fq_concordance_out

    """
        # Split bam file into individual read groups; Ignore MtDNA
        contigs="`samtools view -H input.bam | grep -Po 'SN:([^\\W]+)' | cut -c 4-40 | grep -v 'MtDNA' | tr ' ' '\\n'`"
        rg_list="`samtools view -H input.bam | grep '@RG' | grep -oP 'ID:([^ \\t]+)' | sed 's/ID://g'`"
        samtools split -f '%!.%.' input.bam
        # DO NOT INDEX ORIGINAL BAM; ELIMINATES CACHE!
        bam_list="`ls -1 *.bam | grep -v 'input.bam'`"

        ls -1 *.bam | grep -v 'input.bam' | xargs --verbose -I {} -P ${call_variant_cpus} sh -c "samtools index {}"

        # Call a union set of variants
        for rg in \$rg_list; do
            echo \${contigs} | tr ' ' '\\n' | xargs --verbose -I {} -P ${call_variant_cpus} sh -c "samtools mpileup --redo-BAQ -r {} --BCF --output-tags DP,AD,ADF,ADR,SP --fasta-ref ${reference} \${rg}.bam | bcftools call -T sitelist.tsv.gz --skip-variants indels --multiallelic-caller -O z > {}.\${rg}.vcf.gz"
            order=`echo \${contigs} | tr ' ' '\\n' | awk -v rg=\${rg} '{ print \$1 "." rg ".vcf.gz" }'`
            # Output variant sites
            bcftools concat \${order} -O v | \\
            vk geno het-polarization - | \\
            bcftools query -f '%CHROM\\t%POS[\\t%GT\\t${SM}\\n]' | grep -v '0/1' | awk -v rg=\${rg} '{ print \$0 "\\t" rg }' > \${rg}.rg_gt.tsv
        done;
        cat *.rg_gt.tsv > rg_gt.tsv
        touch out.tsv
        Rscript --vanilla ${fq_concordance_script} 
    """
}

process combine_fq_concordance {

    publishDir analysis_dir + "/concordance", mode: 'copy', overwrite: true

    input:
        file("out*.tsv") from fq_concordance_out.toSortedList()

    output:
        file("fq_concordance.tsv")

    """
        cat <(echo 'a\tb\tconcordant_sites\ttotal_sites\tconcordance\tSM') out*.tsv > fq_concordance.tsv
    """


}

/*
    Call variants
*/

process call_variants_union {

    echo true

    cpus call_variant_cpus

    tag { SM }

    input:
        set SM, file("${SM}.bam"), file("${SM}.bam.bai"), file('sitelist.tsv.gz'), file('sitelist.tsv.gz.tbi') from union_vcf_channel

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
        val merge_vcf from raw_vcf.toSortedList()

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
        set file("RIL.filter.vcf.gz"), file("RIL.filter.vcf.gz.csi") into filtered_vcf

    """

        bcftools view merged.raw.vcf.gz | \\
        vk geno het-polarization - | \\
        bcftools filter --set-GTs . --exclude '((FORMAT/AD[1])/(FORMAT/DP) < 0.75 && FORMAT/GT == "1/1")' - | \\
        bcftools view -O z - > RIL.filter.vcf.gz
        bcftools index -f RIL.filter.vcf.gz
    """
}

filtered_vcf.into { filtered_vcf_stat; hmm_vcf; hmm_vcf_clean; hmm_vcf_out; vcf_tree }


process stat_tsv {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        set file("RIL.filter.vcf.gz"), file("RIL.filter.vcf.gz.csi")  from filtered_vcf_stat

    output:
        file("RIL.filtered.stats.txt")

    """
        bcftools stats --verbose RIL.filter.vcf.gz > RIL.filtered.stats.txt
    """

}

process output_hmm {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        set file("RIL.filter.vcf.gz"), file("RIL.filter.vcf.gz.csi") from hmm_vcf

    output:
        file("gt_hmm.tsv")

    """
        pyenv local anaconda2-4.2.0
        export QT_QPA_PLATFORM=offscreen
        vk hmm --alt=ALT RIL.filter.vcf.gz > gt_hmm.tsv
    """

}

process output_hmm_clean {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        set file("RIL.filter.vcf.gz"), file("RIL.filter.vcf.gz.csi") from hmm_vcf_clean

    output:
        file("gt_hmm_fill.tsv") into gt_hmm_fill
        file("gt_hmm_fill.tsv") into gt_hmm_fill_cross_obj


    """
        pyenv local anaconda2-4.2.0
        export QT_QPA_PLATFORM=offscreen
        vk hmm --infill --endfill --alt=ALT RIL.filter.vcf.gz > gt_hmm_fill.tsv
    """

}



process output_hmm_vcf {

    publishDir analysis_dir + "/vcf", mode: 'copy'

    input:
        set file("RIL.vcf.gz"), file("RIL.vcf.gz.csi") from hmm_vcf_out

    output:
        set file("RIL.hmm.vcf.gz"), file("RIL.hmm.vcf.gz.csi") into gt_hmm

    """
        pyenv local anaconda2-4.2.0
        export QT_QPA_PLATFORM=offscreen
        vk hmm --vcf-out --all-sites --alt=ALT RIL.vcf.gz | bcftools view -O z > RIL.hmm.vcf.gz
        bcftools index RIL.hmm.vcf.gz
    """

}

gt_hmm.into { gt_hmm_tsv; gt_hmm_vcf }

process plot_hmm {

    publishDir analysis_dir + "/hmm", mode: 'copy'

    input:
        file("gt_hmm_fill.tsv") from gt_hmm_fill
        file("script.R") from hmm_plot_script

    output:
        file("gt_hmm.png")
        file("gt_hmm.svg")

    """
        Rscript --vanilla script.R
    """

}

process generate_issue_plots {

    publishDir analysis_dir + "/plots", mode: 'copy'

    input:
        file("SM_bam_idxstats.tsv") from SM_bam_idxstats_plot
        file("fq_coverage.tsv") from fq_coverage_plot
        file("SM_coverage.tsv") from SM_coverage_plot
        file("bam_duplicates.tsv") from bam_duplicates_plot
        file("generate_plots.R") from Channel.fromPath("generate_plots.R")

    output:
        file("coverage_comparison.png")
        file("coverage_comparison.svg")
        file("unmapped_reads.png")
        file("unmapped_reads.svg")
        file("duplicates.png")
        file("duplicates.svg")

    """
        Rscript --vanilla generate_plots.R
    """
}


process output_tsv {

    publishDir analysis_dir + "/hmm", mode: 'copy'

    input:
        set file("RIL.hmm.vcf.gz"), file("RIL.hmm.vcf.gz.csi") from gt_hmm_tsv

    output:
        file("gt_hmm.tsv")

    """
        cat <(echo -e "CHROM\tPOS\tSAMPLE\tGT\tGT_ORIG\tAD") <(bcftools query -f '[%CHROM\t%POS\t%SAMPLE\t%GT\t%GT_ORIG\t%AD\n]' RIL.hmm.vcf.gz | sed 's/0\\/0/0/g' | sed 's/1\\/1/1/g') > gt_hmm.tsv
    """

}

process generate_cross_object {

    publishDir analysis_dir + "/cross_object", mode: 'copy'

    input:
        set file("RIL.hmm.vcf.gz"), file("RIL.hmm.vcf.gz.csi") from gt_hmm_vcf
        file("SM_coverage.tsv") from SM_coverage_cross_obj
        file("gt_hmm_fill.tsv") from gt_hmm_fill_cross_obj

    output:
        file("cross_obj.Rdata")
        file("cross_obj_geno.tsv")
        file("output_strains.tsv")
        file("breakpoint_sites.tsv.gz")

    """

    # Generate breakpoint sites
    cat <(cut -f 1,2 gt_hmm_fill.tsv) <(cut -f 1,3 gt_hmm_fill.tsv) |\
    grep -v 'chrom' |\
    sort -k1,1 -k2,2n |\
    uniq > breakpoint_sites.tsv
    bgzip breakpoint_sites.tsv -c > breakpoint_sites.tsv.gz && tabix -s1 -b2 -e2 breakpoint_sites.tsv.gz

    # Generate output strains list
    awk  '\$2 > 1 && \$2 != "coverage" { print }'  SM_coverage.tsv  | cut -f 1 | sort > cross_obj_strains.tsv

    paste <(echo -e "strain\\t\\t") <(cat cross_obj_strains.tsv| tr '\n' '\t' | sed 's/\t\$//g') > cross_obj_geno.tsv
    bcftools view -T breakpoint_sites.tsv.gz -m 2 -M 2 RIL.hmm.vcf.gz |\
    bcftools query --samples-file cross_obj_strains.tsv -f '%CHROM\\_%POS\\t%CHROM\\t%POS[\\t%GT]\n' |\
    awk  -v OFS='\t' '''
            {   
                gsub("0/0", "N", \$0);
                gsub("1/1", "C", \$0);
                gsub("./.","", \$0);
                gsub("X","Xchr", \$0);
                \$3;
                print
            }
        ''' - >> cross_obj_geno.tsv

    Rscript ${cross_object_script}

    """

}
