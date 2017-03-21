# ril-nf

Align, call variants, and generate datasets for RIAIL sequence data

### Usage

```
# cd to directory of fastqs
nextflow run main.nf
```

__`fq_ril_sheet.tsv`__ is a file that details all of the sequence data. It is generated using the script `scripts/construct_fq_sheet.sh`. `fq_ril_sheet.tsv`

### Generate CB4856 Sitelist

The `CB4856.20160408.sitelist.tsv.gz` file is probably fine forever, but...

```
bcftools view --samples CB4856,N2 -m 2 -M 2 WI.20160408.filtered.vcf.gz | \
vk filter ALT --min=1 - | \
vk filter REF --min=1 - | \
bcftools view --samples CB4856 - | \
bcftools query --include 'FORMAT/GT == "1/1"' -f '%CHROM\t%POS\t%REF,%ALT\n' | \
bgzip -c > CB4856.20160408.sitelist.tsv.gz && tabix -s 1 -b 2 -e 2 CB4856.20160408.sitelist.tsv.gz
```
