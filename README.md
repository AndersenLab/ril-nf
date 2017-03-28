# ril-nf

Align, call variants, and generate datasets for RIAIL sequence data

### Usage

```
# cd to directory of fastqs
nextflow run main.nf -resume
```

__`fq_ril_sheet.tsv`__ is a file that details all of the sequence data. It is generated using the script `scripts/construct_fq_sheet.sh`. `fq_ril_sheet.tsv`

The file has the following format:

| strain   | fastq_pair_id   | library   | fastq-1-path   | fastq-2-path   |
|:-------|:-----------------------|:------------------|:-------------------------------------------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------|
| QX98   | QX98_GCTACGCTGCTCGAA   | GCTACGCTGCTCGAA   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX98_GCTACGCT-GCTCGAA_L003_R1_001.fq.gz   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX98_GCTACGCT-GCTCGAA_L003_R2_001.fq.gz   |
| QX99   | QX99_AGGCAGAAGCTAATC   | AGGCAGAAGCTAATC   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_AGGCAGAA-GCTAATC_L005_R1_001.fq.gz   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_AGGCAGAA-GCTAATC_L005_R2_001.fq.gz   |
| QX99   | QX99_AGGCAGAAGCTAATC   | AGGCAGAAGCTAATC   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_AGGCAGAA-GCTAATC_L006_R1_001.fq.gz   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_AGGCAGAA-GCTAATC_L006_R2_001.fq.gz   |
| QX99   | QX99_CGAGGCTGTGGCAAT   | CGAGGCTGTGGCAAT   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_CGAGGCTG-TGGCAAT_L003_R1_001.fq.gz   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_CGAGGCTG-TGGCAAT_L003_R2_001.fq.gz   |
| QX99   | QX99_CGAGGCTGTGGCAAT   | CGAGGCTGTGGCAAT   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_CGAGGCTG-TGGCAAT_L004_R1_001.fq.gz   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_CGAGGCTG-TGGCAAT_L004_R2_001.fq.gz   |

RIL fastqs are stored in `/projects/b1059/data/fastq/RIL/dna/processed/`.

* Sequence data is very low coverage by design. No pre-processing takes place because variants will be called at specific positions and low quality/adapter contamination are unlikely to be problematic.
* You may wish to perform adapter trimming if contamination is bad. However, most of these reads will be soft-clipped during alignment.

### Generate CB4856 Sitelist

The `CB4856.20160408.sitelist.tsv.gz` file is probably fine forever. It was generated with the following command:

```
bcftools view --samples CB4856,N2 -m 2 -M 2 WI.20160408.filtered.vcf.gz | \
vk filter ALT --min=1 - | \
vk filter REF --min=1 - | \
bcftools view --samples CB4856 - | \
bcftools query --include 'FORMAT/GT == "1/1"' -f '%CHROM\t%POS\t%REF,%ALT\n' | \
bgzip -c > CB4856.20160408.sitelist.tsv.gz && tabix -s 1 -b 2 -e 2 CB4856.20160408.sitelist.tsv.gz
```
