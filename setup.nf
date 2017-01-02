#!/usr/bin/env nextflow
params.directory = '/projects/b1059/data/fastq/WI/dna/processed/**/'
params.analysis_dir = "/projects/b1059/analysis/WI_concordance"
params.kmer_size = 31
params.sketches = 20000
out_tsv = "out_k${params.kmer_size}_s${params.sketches}.tsv"
println "Running Concordance on " + params.directory
println "Final Output will be ${out_tsv}"

Channel.fromFilePairs(params.directory + '*{1,2}P.fq.gz', flat: true)
        .into { fq_pairs }

process dump_json {

    publishDir "/projects/b1059/workflows/concordance-nf", mode: "copy"

    output:
    file 'fq_data.json' into fq_data

    """
        fq dump > fq_data.json
    """
}

process generate_sets {

    module 'R/3.3.1'

    publishDir "/projects/b1059/workflows/concordance-nf", mode: "copy"

    input:
    file 'fq_data.json' from fq_data

    output:
    file 'strain_set.json' into strain_json
    file 'isotype_set.json' into isotype_json

    '''
    #!/usr/bin/env Rscript --vanilla
    library(dplyr)
    fq <- jsonlite::fromJSON("fq_data.json") %>%
      dplyr::filter(strain_type == "WI") %>%
      dplyr::select(original_strain, strain, isotype, library, filename) %>%
      tidyr::unnest(filename) %>%
      dplyr::filter(grepl("b1059", filename)) %>%
      dplyr::filter(grepl("processed", filename))

    # Generate strain and isotype concordance sets
    fstrains <- lapply(split(fq, fq$strain), function(x) {
      unique(x$filename)
    }) %>% jsonlite::toJSON(.)

    readr::write_lines(fstrains, "strain_set.json")


    # Isotype
    fisotype <- lapply(split(fq, fq$isotype), function(x) {
      unique(x$filename)
    })

    readr::write_lines(fstrains, 'isotype_set.json')

    '''
}
