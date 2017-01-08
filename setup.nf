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

    executor 'local'

    publishDir "/projects/b1059/workflows/concordance-nf", mode: "copy"

    output:
    file 'fq_data.json' into fq_data

    """
        fq query 'strain_type = WI, use = True' > fq_data.json
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
    longest_string <- function(s){return(s[which.max(nchar(s))])}
    lcsbstr <- function(a,b) { 
        matches <- gregexpr("M+", drop(attr(adist(a, b, counts=TRUE), "trafos")))[[1]];
        lengths<- attr(matches, 'match.length')
        which_longest <- which.max(lengths)
        index_longest <- matches[which_longest]
        length_longest <- lengths[which_longest]
        longest_cmn_sbstr  <- substring(longest_string(c(a,b)), index_longest , index_longest + length_longest - 1)
        return(longest_cmn_sbstr ) 
    }

    fq <- jsonlite::fromJSON("fq_data.json") %>%
      dplyr::filter(strain_type == "WI") %>%
      dplyr::filter(grepl("b1059", filename)) %>%
      dplyr::filter(grepl("processed", filename)) %>%
      dplyr::select(original_strain, strain, isotype, library, seq_folder, filename) %>%
      tidyr::unnest(filename) %>%
      dplyr::filter(grepl("1P.fq.gz", filename)) %>%
      dplyr::mutate(basename = basename(filename)) %>%
      dplyr::mutate(ID = paste0(seq_folder, "____", gsub("1P.fq.gz", "", basename))) %>%
      dplyr::mutate(LB = library, SM = strain) %>%
      dplyr::mutate(RG = paste0("@RG\\tID:", ID, "\\tLB:", LB, "\\tSM:", SM))

    # Generate strain and isotype concordance sets
    fq_set = list()
    fstrains <- lapply(split(fq, fq$strain), function(i) {
        lapply(split(i, i$RG), function(x) {
            f1 <- x$filename
            f2 <- gsub("1P.fq.gz", "2P.fq.gz", f1)
            set_id <- lcsbstr(basename(f1), basename(f2))
       c(f1, f2, set_id)
       })
    }) %>% jsonlite::toJSON(.)

    readr::write_lines(fstrains, "strain_set.json")


    #fq %>% dplyr::mutate(SM = isotype) %>%
    #dplyr::mutate(RG = paste0("@RG\tID:", ID, "\tLB:", LB, "\tSM:", SM))
    #
    ## Isotype
    #fisotype <- lapply(split(fq, fq$isotype), function(i) {
    #    lapply(split(i, i$RG), function(x) {
    #   c(x$filename, gsub("_1P.fq.gz", "_2P.fq.gz", x$filename))
    #   })
    #}) %>% jsonlite::toJSON(.)

    #readr::write_lines(fstrains, "isotype_set.json")

    '''
}
