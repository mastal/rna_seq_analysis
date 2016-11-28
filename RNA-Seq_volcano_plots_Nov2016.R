################################################
################################################
# Maria Stalteri
# 28/11/2016
# RNA-Seq_volcano_plots_Nov2016.R
#
# tidied up code for volcano plots
# for RNA-Seq data.
#
################################################
################################################
# Code for coloring and labelling points adapted
# from the blog Getting Genetics Done by 
# Stephen Turner.
#
# Volcano Plots code:
# http://www.gettinggeneticsdone.com/2014/05/r-volcano-plots-to-visualize-rnaseq-microarray.html
#
# Code for fancier plots to prevent labels overlapping:
# http://www.gettinggeneticsdone.com/2016/01/repel-overlapping-text-labels-in-ggplot2.html
#
##################################################
##################################################
# 1. Simple B&W plot, no title, print to a file.
#
#    Making volcano plots using results of DESeq2
#    analysis.

library(DESeq2)
library(RCOlorBrewer)
library(BiocInstaller)
biocLite(c("rafalib", "calibrate")
library(rafalib)
library(calibrate)

# The calibrate package is necessary for labelling
# the points.

# Save plots as png files.
# Other formats (jpg, pdf, tiff) are also supported.

png("KC_pirbright_volcano_plain_no_title.png")
plot(deseq2.res.05.Ordered$log2FoldChange,
  -log10(deseq2.res.05.Ordered$pvalue),
    xlab="log2(FoldChange)", 
      ylab="-log10(pvalue)")

dev.off()
# X11cairo 
#       2 

#####################################################
#####################################################
# 2. Simple B&W plot with plot-title.

png("KC_pirbright_volcano_plain.png")
plot(deseq2.res.05.Ordered$log2FoldChange,
  -log10(deseq2.res.05.Ordered$pvalue),
     main = "Volcano Plot for Culicoides Samples, Pirbright Genome",
       xlab="log2(FoldChange)", 
         ylab="-log10(pvalue)",
           xlim=c(-4, 6))

dev.off()

# X11cairo 
#       2 

#####################################################
#####################################################
# 3. Plot with DE genes with adj.p.value < 0.05
#    and absolute value of log2FoldChange > 1.0
#    as red points for upregulated and green for
#    downregulated.

png("KC_pirbright_volcano_rg_points.png")
plot(deseq2.res.05.Ordered$log2FoldChange,
  -log10(deseq2.res.05.Ordered$pvalue),
     main = "Volcano Plot for Culicoides Samples, Pirbright Genome",
       xlab="log2(FoldChange)", 
         ylab="-log10(pvalue)")

with(subset(deseq2.res.05.Ordered, padj< 0.05 & log2FoldChange > 1.0),
    points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

with(subset(deseq2.res.05.Ordered, padj < 0.05 & log2FoldChange < -1.0),
    points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

dev.off()

# X11cairo 
#       2 

########################################################
########################################################
# 4. Plot with vertical lines at log2FoldChange of +1, -1,
#    and horizontal line at 1.30 (-log10(0.05)),
#    the top 20 DE genes by adj.p.value labelled
#    with Gene ID, and red and green points as for 3. above.

# the calibrate library is needed for labelling the points.
library(calibrate)

png("KC_pirbright_volcano_rg_points_labels.png")

plot(deseq2.res.05.Ordered$log2FoldChange,
  -log10(deseq2.res.05.Ordered$pvalue),
     main = "Volcano Plot for Culicoides Samples, Pirbright Genome",
       xlab="log2(FoldChange)", 
         ylab="-log10(pvalue)",
           xlim=c(-4,6))

abline(h=1.30)
abline(v=c(-1,1))

with(subset(deseq2.res.05.Ordered, padj< 0.05 & log2FoldChange > 1.0),
     points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

with(subset(deseq2.res.05.Ordered, padj < 0.05 & log2FoldChange < -1.0),
     points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

# topgenes was an already existing R object,
# the top 20 DE genes ranked by adj.p.value.

# topgenes was generated with this code;
# topgenes <- head(rownames(deseq2.res.05.Ordered), 20)

with(subset(deseq2.res.05.Ordered, padj < 0.05 & abs(log2FoldChange) > 1.0),
    textxy(log2FoldChange, -log10(pvalue),
      labs=topgenes, cex=.8))

dev.off()

# X11cairo 
#       2 

###########################################################
###########################################################

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
 [1] calibrate_1.7.2            MASS_7.3-45               
 [3] BiocInstaller_1.22.3       rafalib_1.0.0             
 [5] RColorBrewer_1.1-2         devtools_1.12.0           
 [7] DESeq2_1.12.4              SummarizedExperiment_1.2.3
 [9] Biobase_2.32.0             GenomicRanges_1.24.2      
[11] GenomeInfoDb_1.8.3         IRanges_2.6.1             
[13] S4Vectors_0.10.3           BiocGenerics_0.18.0       

loaded via a namespace (and not attached):
 [1] genefilter_1.54.2       locfit_1.5-9.1          splines_3.3.1          
 [4] lattice_0.20-33         colorspace_1.2-6        rtracklayer_1.32.2     
 [7] GenomicFeatures_1.24.5  chron_2.3-47            survival_2.39-5        
[10] XML_3.98-1.4            foreign_0.8-66          withr_1.0.2            
[13] DBI_0.5                 BiocParallel_1.6.6      plyr_1.8.4             
[16] zlibbioc_1.18.0         Biostrings_2.40.2       munsell_0.4.3          
[19] gtable_0.2.0            memoise_1.0.0           latticeExtra_0.6-28    
[22] geneplotter_1.50.0      biomaRt_2.28.0          AnnotationDbi_1.34.4   
[25] Rcpp_0.12.6             acepack_1.3-3.3         xtable_1.8-2           
[28] scales_0.4.0            Hmisc_3.17-4            annotate_1.50.1        
[31] XVector_0.12.1          Rsamtools_1.24.0        gridExtra_2.2.1        
[34] ggplot2_2.1.0           digest_0.6.10           grid_3.3.1             
[37] tools_3.3.1             bitops_1.0-6            RCurl_1.95-4.8         
[40] RSQLite_1.0.0           Formula_1.2-1           cluster_2.0.4          
[43] Matrix_1.2-6            data.table_1.9.6        rpart_4.1-10           
[46] GenomicAlignments_1.8.4 nnet_7.3-12            

##########################################################
##########################################################
