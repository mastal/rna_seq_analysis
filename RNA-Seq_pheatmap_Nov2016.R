#############################################################
#############################################################
# Maria Stalteri
# 26/11/2016
# RNA-Seq_pheatmap_nov2016.R
#
# cleaned up sample code for heatmap of top genes by p.adj
# with pheatmap.

#############################################################
#############################################################
# 2. Heatmap using pheatmap function from pheatmap package.
#
#    Code adapted from the RNA-Seq section of the Harvard-X
#    PH525.7x edX course by Rafael Irizarry and Michael Love.
#
#    Workflow starts with the cnts object from the RNA-Seq
#    Workflow in the vignette for the gage Bioconductor
#    package for pathway analysis.
#
#    Use the cnts object after removing rows with zero counts,
#    but before normalization.

library(DESeq2)

# if pheatmap package isn't installed
library(BiocInstaller)
biocLite("pheatmap")
library(pheatmap)

# check the pre-existing counts matrices
dim(hnrnp.cnts)
dim(cnts)
dim(cnts.norm)

# there are also some pre-existing objects used
# with DESeq

coldat
# DataFrame with 6 rows and 1 column
#       grp
#  <factor>
# 1  control
# 2  control
# 3  control
# 4 infected
# 5 infected
# 6 infected

grp.idx
# [1] "control"  "control"  "control"  "infected" "infected" "infected"

# do the DESeq analysis in steps,
# not all at once
dds.step1 <-DESeqDataSetFromMatrix(cnts,
   colData=coldat, design = ~ grp)

dds.step2<-estimateSizeFactors(dds.step1)

# do rlog transformation.
rld <- rlog(dds.step2, blind=FALSE)

rld
# class: DESeqTransform
# ...

###########################################
###########################################
# make heatmap with the top 20 DE genes,
# by adj.p.value
# the deseq2.res.05.ord object had previously been
# created as part of the DESeq2 workflow.

# grp.idx <- rep(c("control", "infected"), each=3)
# coldat=DataFrame(grp=factor(grp.idx))
# dds <- DESeqDataSetFromMatrix(cnts,
#   colData = coldat, design = ~ grp)     
# dds <- DESeq(dds)
# deseq2.res.05 <- results(dds, alpha=0.05)
# deseq2.res.05.ord <- 
#    deseq2.res.05[order(deseq2.res.05$padj),]

############################################
############################################

topgenes <- head(rownames(deseq2.res.05.ord), 20)
mat.rld <- assay(rld)[topgenes,]

class(mat.rld)
# [1] "matrix"

dim(mat.rld)
# [1] 20  6

# take ratio of rlog transformed counts
# to rowMeans of counts
mat.ratios <- mat.rld - rowMeans(mat.rld)

dim(mat.ratios)
# [1] 20  6

# here you need the dds object from the previous
# DESeq2 workflow.
df.rld <- as.data.frame(colData(dds)[, "grp"])

 df.rld
#   colData(dds)[, "grp"]
# 1               control
# 2               control
# 3               control
# 4              infected
# 5              infected
# 6              infected

# the as.data.frame transformation didn't keep the
# row names,
# this could be causing the problems with the plot
# function below.

####################################################
####################################################
# plot the heatmap

pheatmap(mat.ratios, annotation_col=df.rld)

# Error in check.length("fill") :
#  'gpar' element 'fill' must not be length 0

# this works but doesn't give a color code for
# the samples, control and infected
pheatmap(mat.ratios)

####################################################
####################################################
# try to fix the problem.
# 1. give the samples short names so there will
#    be more room for the actual heatmap on the plot.

colnames.mat<-
     c("BF0A", "BF0B", "BF0C", "BF8A","BF8B", "BF8C")

colnames(mat.ratios) <- colnames.mat

# 2. the df.rld data.frame object hasno row names.
#    try assigning it row names

row.names(df.rld)
# [1] "1" "2" "3" "4" "5" "6"

row.names(df.rld)<-colnames.mat

df.rld
#     colData(dds)[, "grp"]
# BF0A               control
# BF0B               control
# BF0C               control
# BF8A              infected
# BF8B              infected

# try the plot again. 
# it works this time.
pheatmap(mat.ratios, annotation_col=df.rld)

# 3. need to give df.rld a shorter col name

colnames(df.rld)<-"group"

df.rld
#        group
# BF0A  control
# BF0B  control
# BF0C  control
# BF8A infected
# BF8B infected
# BF8C infected

# redo the plot, save to a file
pheatmap(mat.ratios, annotation_col=df.rld,
     filename="hamster_pheatmap_top20_by_adjpval.png")

####################################################
####################################################
# suggestions which didn't solve the problem in this case.

library(grid)
get.gpar()

# there were no NA items in the matrix mat.ratios.

####################################################
####################################################
sessionInfo()

R version 3.3.1 (2016-06-21)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: Scientific Linux 7.1 (Nitrogen)

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
[1] pheatmap_1.0.8             DESeq2_1.12.4
[3] SummarizedExperiment_1.2.3 Biobase_2.32.0
[5] GenomicRanges_1.24.2       GenomeInfoDb_1.8.3
[7] IRanges_2.6.1              S4Vectors_0.10.3
[9] BiocGenerics_0.18.0

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.6          RColorBrewer_1.1-2   plyr_1.8.4
 [4] XVector_0.12.1       tools_3.3.1          bitops_1.0-6
 [7] zlibbioc_1.18.0      rpart_4.1-10         RSQLite_1.0.0
[10] annotate_1.50.1      gtable_0.2.0         lattice_0.20-33
[13] Matrix_1.2-6         DBI_0.5              gridExtra_2.2.1
[16] genefilter_1.54.2    cluster_2.0.4        Biostrings_2.40.2
[19] locfit_1.5-9.1       grid_3.3.1           nnet_7.3-12
[22] data.table_1.9.6     AnnotationDbi_1.34.4 XML_3.98-1.4
[25] survival_2.39-5      BiocParallel_1.6.6   foreign_0.8-66
[28] latticeExtra_0.6-28  Formula_1.2-1        geneplotter_1.50.0
[31] ggplot2_2.1.0        Hmisc_3.17-4         Rsamtools_1.24.0
[34] scales_0.4.0         splines_3.3.1        colorspace_1.2-6
[37] xtable_1.8-2         acepack_1.3-3.3      RCurl_1.95-4.8
[40] munsell_0.4.3        chron_2.3-47

####################################################
####################################################
