# usethis::use_r("filename") to create new R script
# Run FastQC on the input file

install.packages("fastqcr")
library(fastqcr)



fqDir <- system.file("testdata", "fastq")
fastqc(fq.dir = fqDir, # FASTQ files directory
       qc.dir = fqDir, # Results direcory
)
