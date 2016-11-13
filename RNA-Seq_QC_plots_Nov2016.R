################################################
################################################
# Maria Stalteri
# 12/11/2016
# RNA-Seq_QC_plots_Nov2016.R
#
# cleaned up example code for QC plots
# for RNA-Seq data.
# 
################################################
################################################
# 1. MA plots using plotMA function from the
#    DESeq2 package. See the package vignette,
#    Section  .

# Use deseq2.res object from the DESeq2 worflow
# (Section 7.1) from the gage vignette,
# and code for plotting from the DESeq2 vignette
# (Section 1.5.1 MA Plot)

library(DESeq2)

?plotMA

# plotMA(object, alpha, main = "",
#  xlab = "mean of normalized counts", ylim, MLE = FALSE, ...)

# note that plotMA uses a cutoff of alpha=0.10
# as default, to use other values, such as 0.05,
# you have to specify alpha=0.05.

# you can also change ylim, default is c(-2,2),
# and you can change the plot title.

# DE genes are plotted as red points,
# points that are outside the range of the graph 
# are shown as open triangles,
# the x-axis is drawn in red.

# save plot as pdf or png file.
# png works with MS-Word and PowerPoint.

# set ylim according to the max log2 fold-change
# in the dataset.

# for this dataset, max lfc was 5.14.
# make plot with horizontal black lines at
# log2 fold-change = +1 and -1 (fold-change of 2).

png("deseq2.res.05.maplot.alpha05.ylim6.lines2blk.KC.2016nov05.png")
plotMA(deseq2.res.05, alpha=0.05, main="Culicoides samples, 
  0h v 24h, DESeq2, p.adj <= 0.05", ylim=c(-6,6))

abline(h=1, lwd=2)
abline(h=-1, lwd=2)
 
dev.off()
# null device 
#           1 

#####################################################
# another png plot, with blue horizontal lines

png("deseq2.res.05.maplot.alpha05.ylim6.lines2blu.KC.2016nov05.png")
plotMA(deseq2.res.05, alpha=0.05, main="Culicoides samples, 
    0h v 24h, DESeq2, p.adj <= 0.05", ylim=c(-6,6))

abline(h=1, lwd=2, col="blue")
abline(h=-1, lwd=2, col="blue")

dev.off()
# null device 
#          1 

#####################################################
# a plot saved as pdf.
# using the default, with no extra horizontal lines.

pdf("deseq2.res.05.maplot.alpha05.ylim6.nolines.KC.2016nov05.pdf")
plotMA(deseq2.res.05, alpha=0.05, main="Culicoides samples, 
 0h v 24h, DESeq2, p.adj <= 0.05", ylim=c(-6,6))

dev.off()
# null device 
#          1 

#####################################################
#####################################################
