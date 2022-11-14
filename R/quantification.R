#' Indexing and quantification using default setting of salmon.
#' More options will come up in future.
#'
#'
#' @export
#'
quantification <- function(species = NA,
                           release = NA,
                           indexName = "index",
                           fastq = NA,
                           quantOut = "quants"){

    # Obtain .cdna.all.fa.gz and .dna.toplevel.fa.gz from Ensembl database
    cdna <- obtaincDNA(species = "Caenorhabditis Elegans", wantedVersion=release)
    dna <- obtainDNA(species = "Caenorhabditis Elegans", wantedVersion=release)

    # Combine the two gz files to one
    rstudioapi::terminalSend(myTerm, paste0("cdna=", cdnaGZ, "\n"))

    rstudioapi::terminalSend(myTerm, paste0("dna=", dnaGZ, "\n"))

    rstudioapi::terminalSend(myTerm, "cat $cdna $dna > gentrome.fa.gz\n")

    # Build index
    rstudioapi::terminalSend(myTerm, paste0("indexName=", indexName, '\n'))

    rstudioapi::terminalSend(myTerm,
                             "salmon index -t gentrome.fa.gz -i $indexName\n")


    # Quantification
    rstudioapi::terminalSend(myTerm, paste0("fastq=", "hello", "\n"))

    rstudioapi::terminalSend(myTerm, paste0("$quantOut=", quantOut, "\n"))

    rstudioapi::terminalSend(myTerm, "salmon quant -i $indexName \
                                      -l A -r $fastq --validateMappings \
                                      -o $quantOut")
}
