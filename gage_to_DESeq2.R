#############################################
#############################################
# Maria Stalteri
# 12/11/2016
# gage_to_DESeq2.R
#
# Cleaned up example of code to take counts table
# from gage and get DE genes with DESeq2.
#
#############################################
#############################################
# This is adapted from the gage RNA-Seq 
# workflow vignette and from the DESeq2
# vignette.
#
# Start with the cnts matrix, after removing
# rows with zero counts with sel.rn
#
##############################################
##############################################
# Start with the gage workflow vignette, Section
# 7.1, DESeq2.

library(DESeq2)

# I think by default the levels get listed in
# alphabetical order, which works out OK here.
grp.idx <-rep(c("control", "infected"), each=3)
 
grp.idx
# [1] "control"  "control"  "control"  "infected" "infected" "infected"

coldat=DataFrame(grp=factor(grp.idx))

# here cnts is the counts matrix from the gage
# workflow, BEFORE normalization, but AFTER
# removing rows with zero counts with [sel.rn, ]
dds <- DESeqDataSetFromMatrix(cnts, 
    colData=coldat, design = ~ grp)

dds<-DESeq(dds)

# note that by default, results() uses
# a cutoff of padj = 0.10, if you want
# to use 0.05, you have to specify it
# with alpha=0.05.
# OUTPUT #1.
deseq2.res.05<-results(dds, alpha=0.05)

# this comes out ordered by the row.names,
# usually the gene IDs, and gives results
# for all genes, not just the DE ones.

deseq2.res.05
# log2 fold change (MAP): grp infected vs control 
# Wald test p-value: grp infected vs control 
# DataFrame with 20338 rows and 6 columns

# write to a text file.


#################################################
#################################################
# this code is taken from the DESeq2 vignette.
summary(deseq2.res.05)

# out of 20338 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 2078, 10% 
# LFC < 0 (down)   : 2024, 10% 
# outliers [1]     : 4, 0.02% 
# low counts [2]   : 6307, 31% 
# (mean count < 17)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

#####################################################
# to get the results ordered by adj.p.value.
# this gives all the genes, not just the ones
# called DE.
# OUTPUT #2.
deseq2.res.05.ord<-
     deseq2.res.05[order(deseq2.res.05$padj),]

deseq2.res.05.ord

# log2 fold change (MAP): grp infected vs control 
# Wald test p-value: grp infected vs control 
# DataFrame with 20338 rows and 6 columns

# write the results to a tab-delimited text file.
write.table(deseq2.res.05.ord, 
  file="Crigri_RS2014_gff3_Entrez_modified_deseq2_05_ord_padj_nov10.txt",
    sep="\t", 
      row.names=TRUE, 
        col.names=NA,
          quote=FALSE)

#################################################
# get a list of only the DE genes, with subset()
# OUTPUT #3.
deseq2.res.05.sig<-
       subset(deseq2.res.05.ord, padj < 0.05)

deseq2.res.05.sig

# log2 fold change (MAP): grp infected vs control 
# Wald test p-value: grp infected vs control 
# DataFrame with 4102 rows and 6 columns

# write results to a text file.
write.table(deseq2.res.05.sig, 
  file="Crigri_RS2014_gff3_Entrez_modified_deseq2_05_sig_only_by_padj_nov10.txt",
    sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

######################################################
# get lists of only upregulated or downregulated,
# ordered by log2 fold-change.
# for the list of upregulated genes,
# you need the negative sign to get the order by
# decreasing lfc, i.e. the largest fold-change first.
# OUTPUT #4.

# this is still ordered by padj
deseq2.res.05.sig.up<-
         subset(deseq2.res.05.sig, log2FoldChange > 0)

# to order by log2 fold-change
deseq2.res.05.sig.up.by.lfc<-
   deseq2.res.05.sig.up[order(-deseq2.res.05.sig.up$log2FoldChange),]

deseq2.res.05.sig.up.by.lfc
# log2 fold change (MAP): grp infected vs control 
# Wald test p-value: grp infected vs control 
# DataFrame with 2078 rows and 6 columns

# write to a file:
write.table(deseq2.res.05.sig.up.by.lfc, 
  file="Crigri_RS2014_gff3_Entrez_modified_deseq2_05_up_by_lfc_nov10.txt",
    sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

#######################################################
# OUTPUT #5.
# just the down-regulated genes.

deseq2.res.05.sig.down
    <-subset(deseq2.res.05.sig, log2FoldChange < 0)

deseq2.res.05.sig.down.by.lfc<-
  deseq2.res.05.sig.down[order(deseq2.res.05.sig.down$log2FoldChange),]

# log2 fold change (MAP): grp infected vs control 
# Wald test p-value: grp infected vs control 
# DataFrame with 2024 rows and 6 columns

# write to a file:

write.table(deseq2.res.05.sig.down.by.lfc, 
  file="Crigri_RS2014_gff3_Entrez_modified_deseq2_05_down_by_lfc_nov10.txt",
    sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

#######################################################
#######################################################
# save info about which package versions were used
# with sessionInfo()

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
[1] DESeq2_1.12.4              SummarizedExperiment_1.2.3
[3] Biobase_2.32.0             GenomicRanges_1.24.2      
[5] GenomeInfoDb_1.8.3         IRanges_2.6.1             
[7] S4Vectors_0.10.3           BiocGenerics_0.18.0       

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.6             RColorBrewer_1.1-2      plyr_1.8.4             
 [4] XVector_0.12.1          GenomicFeatures_1.24.5  bitops_1.0-6           
 [7] tools_3.3.1             zlibbioc_1.18.0         rpart_4.1-10           
[10] biomaRt_2.28.0          annotate_1.50.1         RSQLite_1.0.0          
[13] gtable_0.2.0            lattice_0.20-33         Matrix_1.2-6           
[16] DBI_0.5                 gridExtra_2.2.1         genefilter_1.54.2      
[19] rtracklayer_1.32.2      cluster_2.0.4           Biostrings_2.40.2      
[22] locfit_1.5-9.1          nnet_7.3-12             grid_3.3.1             
[25] data.table_1.9.6        AnnotationDbi_1.34.4    XML_3.98-1.4           
[28] survival_2.39-5         BiocParallel_1.6.6      foreign_0.8-66         
[31] latticeExtra_0.6-28     Formula_1.2-1           geneplotter_1.50.0     
[34] ggplot2_2.1.0           Hmisc_3.17-4            Rsamtools_1.24.0       
[37] scales_0.4.0            GenomicAlignments_1.8.4 splines_3.3.1          
[40] xtable_1.8-2            colorspace_1.2-6        acepack_1.3-3.3        
[43] RCurl_1.95-4.8          munsell_0.4.3           chron_2.3-47         

#######################################################
#######################################################
# save R commands and data objects.

savehistory("pathview_hamster_RS2014_gff3_Entrez_nov10_DESeq2b.Rhistory")
save.image("pathview_hamster_RS2014_gff3_Entrez_nov10_DESeq2b.RData")

#######################################################
#######################################################

