
#Processing FASTQ Files with ShortRead

#install libraries

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("ShortRead")


#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("systemPipeR")


library(ShortRead)
library(systemPipeR)

#place the FASTQ files in a folder named data in WD

fastq <- list.files("data", "*.fastq$"); fastq <- paste("data/", fastq, sep="")
names(fastq) <- paste("flowcell6_lane", 1:length(fastq), sep="_") 
(fq <- readFastq(fastq[1])) # Imports first FASTQ file


countLines(dirPath="./data", pattern=".fastq$")/4 # Counts numbers of reads in FASTQ files


id(fq)[1] # Returns ID field


sread(fq)[1] # Returns sequence


quality(fq)[1] # Returns Phred scores 


as(quality(fq), "matrix")[1:4,1:12] # Coerces Phred scores to numeric matrix


ShortReadQ(sread=sread(fq), quality=quality(fq), id=id(fq)) # Constructs a ShortReadQ from components


library(systemPipeR)
fqlist <- seeFastq(fastq=fastq, batchsize=800, klength=8) # For real data set batchsize to at least 10^5 
seeFastqPlot(fqlist)


