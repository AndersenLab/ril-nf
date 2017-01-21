#!/usr/bin/env nextflow

println params.type

process dump_json {

    executor 'local'

    publishDir "."

    output:
        file "fq_data.json" into fq_data

    """
        fq query 'strain_type = ${params.type}, use = True' > fq_data.json
    """
}

process generate_sets {

    // Remember to add module for R module 'R/3.3.1'

    publishDir ".", mode: 'copy'

    input:
        file "fq_data.json" from fq_data

    output:
        file "strain_set.json" into strain_json

    """
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
      dplyr::filter(strain_type == "${params.type}") %>%
      dplyr::filter(grepl("b1059", filename)) %>%
      dplyr::filter(grepl("processed", filename)) %>%
      dplyr::mutate(library = barcode) %>%
      dplyr::select(strain, library, seq_folder, filename) %>%
      tidyr::unnest(filename) %>%
      dplyr::filter(grepl("1P.fq.gz", filename)) %>%
      dplyr::mutate(basename = basename(filename)) %>%
      dplyr::mutate(ID = paste0(seq_folder, "____", gsub("1P.fq.gz", "", basename))) %>%
      dplyr::mutate(LB = library, SM = strain) %>%
      dplyr::mutate(RG = paste0("@RG\\tID:", ID, "\\tLB:", LB, "\\tSM:", SM))

    # Generate strain and isotype concordance sets
    fstrains <- lapply(split(fq, fq\$strain), function(i) {
        lapply(split(i, i\$RG), function(x) {
            f1 <- x\$filename
            f2 <- gsub("1P.fq.gz", "2P.fq.gz", f1)
            set_id <- lcsbstr(basename(f1), basename(f2))
       c(f1, f2, set_id)
       })
    }) %>% jsonlite::toJSON(., pretty = TRUE)

    readr::write_lines(fstrains, "strain_set.json")
    """
}
