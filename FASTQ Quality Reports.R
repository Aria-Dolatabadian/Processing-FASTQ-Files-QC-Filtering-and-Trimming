
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



#Filtering and Trimming FASTQ Files with ShortRead

fqtrim <- trimLRPatterns(Rpattern="GCCCGGGTAA", subject=fq)
sread(fq)[1:2] # Before trimming


sread(fqtrim)[1:2] # After trimming


#Read counting and duplicate removal 

tables(fq)$distribution # Counts read occurences


sum(srduplicated(fq)) # Identifies duplicated reads


fq[!srduplicated(fq)]


#Trimming low quality tails

cutoff <- 30
cutoff <- rawToChar(as.raw(cutoff+33))
sread(trimTails(fq, k=2, a=cutoff, successive=FALSE))[1:2]


#Removal of reads with Phred scores below a threshold value

cutoff <- 30
qcount <- rowSums(as(quality(fq), "matrix") <= 20) 
fq[qcount == 0] # Number of reads where all Phred scores >= 20


#Removal of reads with x Ns and/or low complexity segments

filter1 <- nFilter(threshold=1) # Keeps only reads without Ns
filter2 <- polynFilter(threshold=20, nuc=c("A","T","G","C")) # Removes reads with nucleotide bias, >=20 of any base
filter <- compose(filter1, filter2)
fq[filter(fq)]



#Memory Efficient FASTQ Processing
#Streaming through FASTQ files with FastqStreamer and random sampling reads with FastqSampler

fq <- yield(FastqStreamer(fastq[1], 50)) # Imports first 50 reads 
fq <- yield(FastqSampler(fastq[1], 50)) # Random samples 50 reads 


#Streaming through a FASTQ file while applying filtering/trimming functions and writing the results to a new file here SRR038845.fastq_sub in data directory.

f <- FastqStreamer(fastq[1], 50) 
while(length(fq <- yield(f))) {
    fqsub <- fq[grepl("^TT", sread(fq))] 
    writeFastq(fqsub, paste(fastq[1], "sub", sep="_"), mode="a", compress=FALSE)
}
close(f)





