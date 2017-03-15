#!/usr/bin/bash
# Generate fq_sheet for 150225_D00422_0165_BHGMN3ADXX
prefix="/projects/b1059/data/fastq/RIL/dna/processed/150225_D00422_0165_BHGMN3ADXX"
cat ${prefix}/.fqdata | grep "1P" | awk -v prefix=${prefix} '{    
                                    fq2 = $1; ID = $1; LB = $4;
                                    gsub("1P.fq.gz", "2P.fq.gz", fq2);
                                    gsub("_1P.fq.gz", "", ID);
                                    print $2 "\t" ID "\t" LB "\t" prefix "/" $1 "\t" prefix "/" fq2;
                                }' >> fq_temp.tsv

# Generate fq_sheet for 150721_700819F_0366_AHTY7WADXX-ECA1-4
prefix="/projects/b1059/data/fastq/RIL/dna/processed/150721_700819F_0366_AHTY7WADXX-ECA1-4"
ls -1 ${prefix} | grep '_R1_' | awk -v prefix=${prefix} '{
                            split($0, a, "_");
                            SM = a[1]; 
                            gsub("-","",a[2]); // LB
                            ID = a[1] "_" a[2];
                            fq2 = $1; gsub("_R1_", "_R2_", fq2);
                            print SM "\t" ID "\t" a[2] "\t" prefix "/" $1 "\t" prefix "/" fq2;
                            }' >> fq_temp.tsv


# Generate fq_sheet for 151009_D00422_0262_BC7NJ0ANXX-ECA
prefix="/projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA"
ls -1 ${prefix} | grep '_R1_' | awk -v prefix=${prefix} '{
                            split($0, a, "_");
                            SM = a[1]; 
                            gsub("-","",a[2]); // LB
                            ID = a[1] "_" a[2];
                            fq2 = $1; gsub("_R1_", "_R2_", fq2);
                            print SM "\t" ID "\t" a[2] "\t" prefix "/" $1 "\t" prefix "/" fq2;
                            }' >> fq_temp.tsv


cat fq_temp.tsv | sort -k1,1 > ../fq_ril_sheet.tsv
rm fq_temp.tsv