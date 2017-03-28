# ril-nf

Align, call variants, and generate datasets for RIAIL sequence data

## Usage

```
# cd to directory of fastqs
nextflow run main.nf -resume
```

__`fq_ril_sheet.tsv`__ is a file that details all of the sequence data. It is generated using the script `scripts/construct_fq_sheet.sh`. RIL fastqs are stored in `/projects/b1059/data/fastq/RIL/dna/processed/`. The script `scripts/construct_fq_sheet.sh` uses the directory structure and filenames to construct the `fq_ril_sheet.tsv` file.

The file has the following format:

| strain   | fastq_pair_id   | library   | fastq-1-path   | fastq-2-path   |
|:-------|:-----------------------|:------------------|:-------------------------------------------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------|
| QX98   | QX98_GCTACGCTGCTCGAA   | GCTACGCTGCTCGAA   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX98_GCTACGCT-GCTCGAA_L003_R1_001.fq.gz   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX98_GCTACGCT-GCTCGAA_L003_R2_001.fq.gz   |
| QX99   | QX99_AGGCAGAAGCTAATC   | AGGCAGAAGCTAATC   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_AGGCAGAA-GCTAATC_L005_R1_001.fq.gz   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_AGGCAGAA-GCTAATC_L005_R2_001.fq.gz   |
| QX99   | QX99_AGGCAGAAGCTAATC   | AGGCAGAAGCTAATC   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_AGGCAGAA-GCTAATC_L006_R1_001.fq.gz   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_AGGCAGAA-GCTAATC_L006_R2_001.fq.gz   |
| QX99   | QX99_CGAGGCTGTGGCAAT   | CGAGGCTGTGGCAAT   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_CGAGGCTG-TGGCAAT_L003_R1_001.fq.gz   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_CGAGGCTG-TGGCAAT_L003_R2_001.fq.gz   |
| QX99   | QX99_CGAGGCTGTGGCAAT   | CGAGGCTGTGGCAAT   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_CGAGGCTG-TGGCAAT_L004_R1_001.fq.gz   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_CGAGGCTG-TGGCAAT_L004_R2_001.fq.gz   |


* Sequence data is very low coverage by design. No pre-processing takes place because variants will be called at specific positions and low quality/adapter contamination are unlikely to be problematic.
* You may wish to perform adapter trimming if contamination is bad. However, most of these reads will be soft-clipped during alignment.

## Generate CB4856 Sitelist

The `CB4856.20160408.sitelist.tsv.gz` file is probably fine forever. It was generated with the following command:

```
bcftools view --samples CB4856,N2 -m 2 -M 2 WI.20160408.filtered.vcf.gz | \
vk filter ALT --min=1 - | \
vk filter REF --min=1 - | \
bcftools view --samples CB4856 - | \
bcftools query --include 'FORMAT/GT == "1/1"' -f '%CHROM\t%POS\t%REF,%ALT\n' | \
bgzip -c > CB4856.20160408.sitelist.tsv.gz && tabix -s 1 -b 2 -e 2 CB4856.20160408.sitelist.tsv.gz
```

## Output

The output directory looks like this:

```
├── concordance
│   └── fq_concordance.tsv
├── cross_object
│   ├── breakpoint_sites.tsv.gz
│   ├── cross_obj.Rdata
│   ├── cross_obj_geno.tsv
│   ├── cross_obj_pheno.tsv
│   └── cross_obj_strains.tsv
├── duplicates
│   └── bam_duplicates.tsv
├── fq
│   ├── fq_bam_idxstats.tsv
│   ├── fq_bam_stats.tsv
│   ├── fq_coverage.full.tsv
│   └── fq_coverage.tsv
├── fq_ril_sheet.tsv
└── hmm
    ├── gt_hmm.png
    ├── gt_hmm.svg
    └── gt_hmm.tsv
```

### concordance/
 
__fq_concordance.tsv__ - Contains the concordance among fastqs belonging to the same strain. 

* `a` - First fastq
* `b` - Second fastq
* `concordant_sites` - Number of sites concordant between `a` and `b`
* `total_sites` - The total number of sites called in both `a` and `b`
* `concordance` - frequency concordant.
* `SM` - Strain

| a                     | b                     |   concordant_sites |   total_sites |   concordance | SM    |
|:----------------------|:----------------------|-------------------:|--------------:|--------------:|:------|
| QX204_CAGAGAGGCCGGATA | EA-G02_TTAACTC_L001   |             148794 |        151385 |      0.982885 | QX204 |
| EA-G02_TTAACTC_L001   | QX204_CAGAGAGGCCGGATA |             148794 |        151385 |      0.982885 | QX204 |

### cross_object/

* __breakpoint_sites.tsv.gz__ - A gzipped file of sites flanking loci where recombination has occured.
* __cross_obj.Rdata__ - A cross object usable with `qtl`
* __cross_obj_geno.tsv__ - Genotypes file used by the R `qtl` package.
* __cross_obj_pheno.tsv__ - Phenotypes file used by the R `qtl` package.
* __cross_obj_strains.tsv__ - List of strains within the cross object.

The script used to generate the cross object is called `generate_cross_object.R` and is located in the root of this repo.

### duplicates/

__bam_duplicates.tsv__ - A summary of duplicate reads from aligned bams.

### fq/

* __fq_bam_idxstats.tsv__
* __fq_bam_stats.tsv__
* __fq_coverage.full.tsv__
* __fq_coverage.tsv__

### SM/

* __SM_bam_idxstats.tsv__
* __SM_bam_stats.tsv__
* __SM_coverage.full.tsv__
* __SM_coverage.tsv__

### hmm/

* __gt_hmm.(png/svg)__ - Image depicting RILs.
* __gt_hmm.tsv__ - Long form genotypes file.

### plots/

__coverage_comparison.(png/svg)__

__duplicates.(png/svg)__

__unmapped_reads.(png/svg)__

### vcf/

* __gt_hmm.tsv__ - Haplotypes defined by region with associated information. 
* __gt_hmm_fill.tsv__ - Same as above, but using `--infill` and `--endfill` with VCF-Kit. For more information, see [VCF-Kit DOcumentation](http://vcf-kit.readthedocs.io/en/latest/)
* __RIL.filter.vcf.gz__ - A VCF of filtered genotypes. 
* __RIL.filtered.stats.txt__ - Summary of filtered genotypes. Generated by `bcftools stats RIL.filter.vcf.gz`
* __RIL.hmm.vcf.gz__ - The RIL VCF as output by VCF-Kit; HMM applied to determine genotypes.
* __union_vcfs.txt__ - A list of VCFs that were merged to generate RIL.filter.vcf.gz
