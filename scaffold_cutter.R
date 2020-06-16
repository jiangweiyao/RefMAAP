#!/usr/bin/env Rscript
## Usage is scaffold_cutter.R <topN> <fasta_file> <depth_file> <out_fasta> <out_table>

library(Biostrings)

args = commandArgs(trailingOnly=TRUE)

topN = as.numeric(args[1])
fasta_file = args[2]
depth_file = args[3]
out_fasta = args[4]
out_table = args[5]

depth_count <- read.table(file = depth_file, sep = "\t")
colnames(depth_count) <- c("RefSeq", "Position", "Depth")
depth_count$Range <- (depth_count$Position + 1 == c(depth_count$Position[-1], 0)) & c(FALSE, head((depth_count$Position + 1 == c(depth_count$Position[-1], 0)), -1))
depth_count_boundary <- subset(depth_count, Range == FALSE)

range_list  <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(range_list) <- c("RefSeq", "Begin", "End", "AveDepth")

for(i in seq(1, dim(depth_count_boundary)[1], 2)){
    RefSeqi=as.character(depth_count_boundary[i,"RefSeq"])
    Begini=as.numeric(depth_count_boundary[i,"Position"])
    Endi=as.numeric(depth_count_boundary[i+1,"Position"])
    AveDepthi=mean(subset(depth_count, (RefSeq == RefSeqi)&(Position >= Begini)& (Position <= Endi))$Depth)
    range_summary <- data.frame(RefSeq=RefSeqi,
                                Begin=Begini,
                                End=Endi,
                                AveDepth=AveDepthi)
    range_list <- rbind(range_list, range_summary)
}

RefSeq_coverage <- aggregate(range_list$AveDepth, by=list(range_list$RefSeq), FUN=mean)
RefSeq_coverage <- RefSeq_coverage[order(order(RefSeq_coverage[,2], decreasing = TRUE)),]

topN_ref <- as.character(head(RefSeq_coverage, topN)[,1])
range_list_topN <- subset(range_list, RefSeq %in% topN_ref)
write.csv(range_list_topN, out_table, row.names = FALSE)
   
fasta <- readDNAStringSet(fasta_file)
fasta_list <- DNAStringSet()

for(i in seq(1, dim(range_list_topN)[1], 1)){
    fastai <- DNAStringSet(fasta[[as.character(range_list_topN[i, 1])]][range_list_topN[i, 2]:range_list_topN[i, 3]])
    #print(i)
    fasta_list <- c(fasta_list, fastai)
    #print(fasta_list)
    #writeXStringSet(fastai, "fasta.fasta", append=TRUE)
}

names(fasta_list) <- paste(as.character(range_list_topN[, 1]), as.character(range_list_topN[, 2]), as.character(range_list_topN[, 3]), sep = "_")
writeXStringSet(fasta_list, out_fasta)
