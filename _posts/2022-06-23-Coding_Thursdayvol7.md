---
title: TT-seq analysis - Scale factors
date: 2022-06-23 
categories: ["Coding Thursday"]
background: /assets/theme/images/ttseq.png
---

In this post, I will normalize the sequencing samples, using the read counts from ERCC spike-ins. To this end I will generate scale factors, with the assumption that spike-ins are equally represented in each sample. Spike-in counts are stored in `*.ReadsPerGene.out.tab` files generated during samples alignment to spike-in sequences. 

1.Creating count matrices

```bash
library(dplyr)
dirPath <- "/your/folder/with/*.ReadsPerGene.out.tab"
setwd(dirPath)
# only return file names with a given pattern
dir(pattern="ReadsPerGene.out.tab")
## [1] "SRR13866843_GSM5137573_FRNA_2i_2d_rep1.ReadsPerGene.out.tab"
## [2] "SRR13866844_GSM5137574_FRNA_2i_2d_rep2.ReadsPerGene.out.tab"
## [3] "SRR13866845_GSM5137575_FRNA_2i_7d_rep1.ReadsPerGene.out.tab"
## [4] "SRR13866853_GSM5137583_LRNA_2i_2d_rep1.ReadsPerGene.out.tab"
## [5] "SRR13866854_GSM5137584_LRNA_2i_2d_rep2.ReadsPerGene.out.tab"
## [6] "SRR13866855_GSM5137585_LRNA_2i_7d_rep1.ReadsPerGene.out.tab"

# save the results to a variable
files <- dir(pattern="ReadsPerGene.out.tab")

#create count matrix for all samples, rows are spikes, columns are samples
counts <- c()
for( i in seq_along(files) ){
  x <- read.table(file=files[i], sep="\t", header=F, as.is=T)
  counts <- cbind(counts, x[,2])
}

# set the row names
rownames(counts) <- x[,1]
# set the column names based on input file names, with pattern removed
colnames(counts) <- sub("_ReadsPerGene.out.tab","",files)
#select rows with specific rownames
# specify rows to keep
# based on Shao et al. total transcriptome (F) has 6 spike ins, whilst labeled (L) has only 4
keep_in_F<-c("ERCC-00043","ERCC-00136","ERCC-00145","ERCC-00092","ERCC-00002","ERCC-00170")
keep_in_L<-c("ERCC-00043","ERCC-00136","ERCC-00145","ERCC-00092")

#seperate the total from labelled samples, FRNA, LRNA
spike_matrix_total<-(counts[rownames(counts) %in% keep_in_F,])[,grepl( "FRNA" , colnames(counts))]
spike_matrix_label<-(counts[rownames(counts) %in% keep_in_L,])[,grepl( "LRNA" , colnames(counts))]
#change column names for readability
colnames(spike_matrix_label)<-c("LRNA_2i_2d_rep1","LRNA_2i_2d_rep2", "LRNA_2i_7d_rep1")
colnames(spike_matrix_total)<-c("FRNA_2i_2d_rep1","FRNA_2i_2d_rep2","FRNA_2i_7d_rep1")
write.table(spike_matrix_label, file="spikematrix_L.txt", row.names=TRUE, col.names=TRUE)
write.table(spike_matrix_total, file="spikematrix_F.txt", row.names=TRUE, col.names=TRUE)
```

2. Calculating scale factor using estimateSizeFactors() function in the DESeq2 package (Love et al.)

```{r}
library(DESeq2)
## Loading required package: S4Vectors
## Loading required package: stats4
## Loading required package: BiocGenerics
## Warning: package 'BiocGenerics' was built under R version 4.0.5
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## The following objects are masked from 'package:dplyr':
## 
##     combine, intersect, setdiff, union
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which.max, which.min
## 
## Attaching package: 'S4Vectors'
## The following objects are masked from 'package:dplyr':
## 
##     first, rename
## The following object is masked from 'package:base':
## 
##     expand.grid
## Loading required package: IRanges
## 
## Attaching package: 'IRanges'
## The following objects are masked from 'package:dplyr':
## 
##     collapse, desc, slice
## Loading required package: GenomicRanges
## Loading required package: GenomeInfoDb
## Warning: package 'GenomeInfoDb' was built under R version 4.0.5
## Loading required package: SummarizedExperiment
## Loading required package: MatrixGenerics
## Loading required package: matrixStats
## 
## Attaching package: 'matrixStats'
## The following object is masked from 'package:dplyr':
## 
##     count
## 
## Attaching package: 'MatrixGenerics'
## The following objects are masked from 'package:matrixStats':
## 
##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##     colWeightedMeans, colWeightedMedians, colWeightedSds,
##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##     rowWeightedSds, rowWeightedVars
## Loading required package: Biobase
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Attaching package: 'Biobase'
## The following object is masked from 'package:MatrixGenerics':
## 
##     rowMedians
## The following objects are masked from 'package:matrixStats':
## 
##     anyMissing, rowMedians

#DESeq requires data with sample information. Make sure that sample names match the names in count matrix and that the samples are
in the right order. DESeq2 will output an error if this is not the case.

samples_L<-c("LRNA_2i_2d_rep1","LRNA_2i_2d_rep2", "LRNA_2i_7d_rep1")
cond_L<-c("2d","2d","7d")
meta_L<-as.matrix(data.frame(samples_L,cond_L))

#create DEseq object
dds_L<-DESeqDataSetFromMatrix(countData=spike_matrix_label,
colData=meta_L,
design= ~ cond_L)

#Estimate size factors
dds_L <- estimateSizeFactors(dds_L)
sizeFactors(dds_L)

## LRNA_2i_2d_rep1 LRNA_2i_2d_rep2 LRNA_2i_7d_rep1 
##       1.0239437       0.9022947       1.0188438

#repeat for total RNA
samples_F<-c("FRNA_2i_2d_rep1","FRNA_2i_2d_rep2","FRNA_2i_7d_rep1")
cond_F<-c("2d","2d","7d")
meta_F<-as.matrix(data.frame(samples_F,cond_F))
dds_F<-DESeqDataSetFromMatrix(countData=spike_matrix_total,
colData=meta_F,
design= ~ cond_F)

dds_F <- estimateSizeFactors(dds_F)
sizeFactors(dds_F)
## FRNA_2i_2d_rep1 FRNA_2i_2d_rep2 FRNA_2i_7d_rep1 
##       0.7051597       0.5022989       2.8893677
```

3. References

Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.
