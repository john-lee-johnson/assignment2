##Assignment2 is to redo the analysis of Assignment1 using R exclusively
##First step is to read in bed files
##Second step is to find overlapping genes
##Third step is to generate TSS and add 5kbp around each TSS
##Fourth step is to find overlapping marks and TSS
##Fifth step is to find overlapping marks

##Loading Libraries
library(GenomicRanges)
library(rtracklayer)
library(RColorBrewer)
library(compare)
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(IRanges)

#Read bed files
k27ac <- import.bedGraph(con="K27Ac_D7_Th1-W200-G200-E200-island.bedgraph", asRangedData=FALSE)
k4m3 <- import.bedGraph(con="K4m3_Th1_72h-W200-G200-E200-island.bedgraph", asRangedData=FALSE)
refgene <- import.bed(con="mm9_RefSeq.bed", asRangedData=FALSE)
sort1 <- import.bed("1sort.bed")
sort2 <- import.bed("2sort.bed")

#Subsetting the overlapping genes and counting unique genes
k27ac.gene.overlaps <- subsetByOverlaps(refgene, k27ac)
list1 <- sort(unique(k27ac.gene.overlaps$name))
length(intersect(list1,sort1$name))/length(list1)

k4m3.gene.overlaps <- subsetByOverlaps(refgene, k4m3)
list2 <- sort(unique(k4m3.gene.overlaps$name))
length(intersect(list2,sort2$name))/length(list2)

##Creating tss start sites and finding histone mark overlaps
tss <- flank(refgene, width = 5000, both=TRUE, start=TRUE)
k27ac.tss.overlaps <- subsetByOverlaps(tss, k27ac)
k27ac.tss <- subsetByOverlaps(k27ac, tss)
list3 <- sort(unique(k27ac.tss.overlaps$name))
k4m3.tss.overlaps <- subsetByOverlaps(tss, k4m3)
k4m3.tss <- subsetByOverlaps(k4m3, tss)
list4 <- sort(unique(k4m3.tss.overlaps$name))

##Finding overlaps between histone marks
k27ac.overlaps <- subsetByOverlaps(k4m3, k27ac)
k4m3.overlaps <- subsetByOverlaps(k27ac, k4m3)

##Calculating percentages
sum(width(reduce(k27ac.tss, ignore.strand=TRUE)))/sum(width(k27ac))
sum(width(reduce(k4m3.tss, ignore.strand=TRUE)))/sum(width(k4m3))
sum(width(reduce(k27ac.overlaps, ignore.strand=TRUE)))/sum(width(k4m3))
sum(width(reduce(k4m3.overlaps, ignore.strand=TRUE)))/sum(width(k27ac))
