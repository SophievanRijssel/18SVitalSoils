rm(list=ls())
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

setDadaOpt(OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16)
getwd()
project <- "~/18sChrono/"

data <- paste(project, "Raw_data/",sep="")
path <- paste(project, "analysis_18S_PR2/",sep="")
list.files(data)

fnFs <- sort(list.files(data, pattern="_1.fq.gz", full.names = TRUE))
#fnRs <- sort(list.files(path, pattern="_2.fq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), '[', 1)

#library(qckitfastq)


# fseq <- seqTools::fastqq(fnFs[1],k=15)
# read_length(fseq)
# 
#plotComplexity(fnFs[1:20])
#plotQualityProfile(fnFs[1:20])
# plotQualityProfile(fnFs[1:20])
# myPlots<-list()
# pdf("plots18SChronoQualityProfile.pdf",onefile = TRUE)
# for (i in 1:449){
#   myPlots[[i]]<-plotQualityProfile(fnFs[i])
#   print(myPlots[[i]])
# }
# dev.off()

FWD <- "GGCAAGTCTGGTGCCAG"  ## 3NDf forward primer Cavalier-Smith 2009
REV <- "TCCGTCAATTYCTTTAAGT"  ## 1132rmod backward primer Geisen 2018

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names
out<-filterAndTrim(fnFs, filtFs, maxN = 0, multithread = TRUE) #default; truncQ=2; only n's are removed. nothing else done.
head(out)
View(out)
#this part is only to check if all primers are removed
# allOrients <- function(primer) {
#   # Create all orientations of the input sequence
#   require(Biostrings)
#   dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
#   orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
#                RevComp = reverseComplement(dna))
#   return(sapply(orients, toString))  # Convert back to character vector
# }
# FWD.orients <- allOrients(FWD)
# REV.orients <- allOrients(REV)
# FWD.orients
# 
# #fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
# #fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
# # outFreddy.mintrunclen<-filterAndTrim(fnFs, filtFs, maxN=0, maxEE=c(5), truncQ=2, rm.phix=TRUE,
# #                                      compress=TRUE, multithread=20)
# # outFreddy.plustrunclen<-filterAndTrim(fnFs, filtFs, maxN=0, maxEE=c(5), truncQ=2, rm.phix=TRUE,
# #                                      compress=TRUE, multithread=20)
# # outadadaptedMattias<-filterAndTrim(fnFs, filtFs, maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
# #                                       compress=TRUE, multithread=20)
# #truncQ; staat laag zodat ie leert (via expected errors), dat klopt via de manual (dus hier niet trimmen op kwaliteit)
# # head(outITSexample)
# # head(outFreddy.mintrunclen)
# # head(outFreddy.plustrunclen)
# # head(outadadaptedMattias)
# 
# primerHits <- function(primer, fn) {
#   # Counts number of reads in which the primer is found
#   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
#   return(sum(nhits > 0))
# }
# 
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = filtFs[[1]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = filtFs[[1]]))
# 
# 
# FWD.ForwardReads<-as.data.frame(matrix(nrow = 449, ncol = 4))
# REV.ForwardReads<-as.data.frame(matrix(nrow = 449, ncol = 4))
# for (i in 1:449){
#   #FWD.ForwardReads[i,] = sapply(FWD.orients, primerHits, fn = filtFs[[i]])
#   REV.ForwardReads[i,] = sapply(REV.orients, primerHits, fn = filtFs[[i]])
# }
# View(REV.ForwardReads)

#what are the results from above? complementary reverse should be in the raw, but not in the filtered data
#In Freddys 0-51 complementary reverse reads, max 15 if maxEE is set to 2 (more strict than 5), but still not completely fixed; therefore we use cutadapt
#then: select short reads and assign taxonomy..??? 


#https://benjjneb.github.io/dada2/ITS_workflow.html > remove primers by cutadapt
cutadapt <- "/home/nioo/sophier/.conda/envs/cutadaptenv/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version")
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
#fnRs.cut <- file.path(path.cut, basename(fnRs))


FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# # Run Cutadapt
# for(i in seq_along(fnFs)) {
#   system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
#                              "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
#                              fnFs.filtN[i], fnRs.filtN[i])) # input files
# }

#only forward:
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i],  # output files
                             filtFs[i])) # input files
}
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]))

#####vervolg Freddys pipeline#####
# outITS <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
#                      truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
# head(out)
cutFs <- sort(list.files(path.cut, pattern = "_1.fq.gz", full.names = TRUE))
# get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
# sample.names <- unname(sapply(cutFs, get.sample.name))
names(cutFs) <- sample.names
#fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
filtFsmax150 <- file.path(path, "filteredmax150", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFsmax150) <- sample.names
filtFs150240 <- file.path(path, "filtered150240", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs150240) <- sample.names
filtFs240 <- file.path(path, "filtered240", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs240) <- sample.names

outmax150 <- filterAndTrim(cutFs, filtFsmax150, maxN = 0, maxEE = c(2),
                     truncQ = 2, maxLen = 150, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(outmax150)

out150240 <- filterAndTrim(cutFs, filtFs150240, maxN = 0, maxEE = c(2),
                           truncQ = 2, minLen = 150, maxLen = 240, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out150240)

out240 <- filterAndTrim(cutFs, filtFs240, truncLen=c(240),
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=20) # On Windows set multithread=FALSE
head(out240)

# #Freddy
# out <- filterAndTrim(fnFs, filtFs, truncLen=c(240),
#                      maxN=0, maxEE=c(5), truncQ=2, rm.phix=TRUE,
#                      compress=TRUE, multithread=20) # On Windows set multithread=FALSE
# head(out)

#out is afterwards only used in track -> does it do it within the filtFnFs??
datafilt240 <- paste(project, "analysis_18S_PR2/filtered240",sep="")
list.files(datafilt240)
filtFs240 <- sort(list.files(datafilt240, full.names = TRUE))

errF240 <- learnErrors(filtFs240, multithread=20)
plotErrors(errF240, nominalQ=TRUE)

derepFs240 <- derepFastq(filtFs240, verbose = TRUE)
dadaFs240 <- dada(derepFs240, err=errF240, multithread=20)

#dadaFs240 <- dada(filtFs240, err=errF240, multithread=FALSE) #doesnt work

seqtab240 <- makeSequenceTable(dadaFs240)
dim(seqtab240)
table(nchar(getSequences(seqtab240)))

seqtab240.nochim <- removeBimeraDenovo(seqtab240, method="consensus", multithread=20, verbose=TRUE)
dim(seqtab240.nochim)
sum(seqtab240.nochim)/sum(seqtab)


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs240, getN), rowSums(seqtab240.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)

# getN <- function(x) sum(getUniques(x))
# track <- cbind(out, sapply(dadaFs240, getN), rowSums(seqtab240.nochim))
# colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
# rownames(track) <- sample.names
# head(track)
write.csv(track,"overviewreadsduringpipeline.csv")

taxa240 <- assignTaxonomy(seqtab240.nochim, "~/18sChrono/pr2_version_4.14.0_SSU_dada2.fasta.gz", multithread=20)
taxa.print <- taxa240 # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(seqtab240.nochim, "seqtab_sophie-18S_PR2_240.rds")
saveRDS(taxa240, "taxa_sophie-18S_PR2_240.rds")

# /mnt/nfs/bioinfdata/home/NIOO/freddyh/maddy_18S/seqtab_maddy-18S_100.rds


write.csv(seqtab240.nochim, "seqtab_sophie-18S_PR2_240.csv")
write.csv(taxa240, "taxa_sophie-18S_PR2_240.csv")

#####max150#####
# #remove files that do not exist anymore
# datafiltmax150 <- paste(project, "analysis_18S_PR2/filteredmax150",sep="")
# list.files(datafiltmax150)
# filtFsmax150 <- sort(list.files(datafiltmax150, full.names = TRUE))
# 
# errF150 <- learnErrors(filtFsmax150, multithread=20,)
# plotErrors(errF150, nominalQ=TRUE)
# 
# derepFsmax150 <- derepFastq(filtFsmax150, verbose = TRUE)
# dadaFsmax150 <- dada(derepFsmax150, err=errF150, multithread=20)
# 
# #dadaFsmax150 <- dada(filtFsmax150, err=errF150, multithread=20)
# 
# seqtabmax150 <- makeSequenceTable(dadaFsmax150)
# dim(seqtabmax150)
# table(nchar(getSequences(seqtabmax150)))
# 
# seqtabmax150 <- makeSequenceTable(FiltFsmax150)
# dim(seqtabmax150)
# table(nchar(getSequences(seqtabmax150)))
# 
# 
# seqtabmax150.nochim <- removeBimeraDenovo(seqtabmax150, method="consensus", multithread=20, verbose=TRUE)
# dim(seqtabmax150.nochim)
# sum(seqtabmax150.nochim)/sum(seqtabmax150)
# 
# getN <- function(x) sum(getUniques(x))
# track <- cbind(out, sapply(dadaFsmax150, getN), rowSums(seqtabmax150.nochim))
# colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
# rownames(track) <- sample.names
# head(track)
# 
# #https://github.com/pr2database/pr2database/releases #get pr2 database (v 4.14 now!)
# taxamax150 <- assignTaxonomy(seqtabmax150.nochim, "~/18sChrono/pr2_version_4.14.0_SSU_dada2.fasta.gz", multithread=20) 
# taxa.print <- taxamax150 # Removing sequence rownames for display only
# rownames(taxa.print) <- NULL
# head(taxa.print)
# 
# saveRDS(seqtabmax150.nochim, "seqtab_sophie-18S_PR2_max150.rds")
# saveRDS(taxamax150, "taxa_sophie-18S_PR2_max150.rds")
# 
# write.csv(seqtabmax150.nochim, "seqtab_sophie-18S_PR2_max150.csv")
# write.csv(taxamax150, "taxa_sophie-18S_PR2_max150.csv")


