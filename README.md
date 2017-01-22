# concordance-nf

Perform Fastq-Profiling in the current working directory.

### Usage

```
# cd to directory of fastqs
nextflow run Andersenlab/concordance-nf
```

### Generate CB4856 Sitelist
```
bcftools view --samples CB4856 -m 2 -M 2 WI.20160408.filtered.vcf.gz | \
vk filter ALT --min=1 - | \
bcftools query --include 'FORMAT/GT == "1/1"' -f '%CHROM\t%POS\t%REF,%ALT\n' | \
bgzip -c > CB4856.20160408.sitelist.tsv.gz && tabix -s 1 -b 2 -e 2 CB4856.20160408.sitelist.tsv.gz
```
