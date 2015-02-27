###################################################
### code chunk number : library
###################################################
pdf("plots.pdf")
library(S4Vectors)
library(GenomicRanges)

###################################################
### code chunk number : banana experiment design
###################################################
dir <- "/Users/jason/banana/bam"
sampleTable <- readRDS("sample_table_better.RDS")

###################################################
### code chunk number : alignment samfiles
###################################################

# library(Rsamtools)
# filenames <- file.path(dir,paste0(sampleTable$SampleName, ".bam"))
# bamfiles <- BamFileList(filenames, yieldSize=2000000)

###################################################
### code chunk number : sequence features
###################################################
# library(GenomicFeatures)
# gtffile <- file.path(dir, "../assembly","musa.gtf")
# txdb <- makeTranscriptDbFromGFF(gtffile, format = "gtf")
# genes <- exonsBy(txdb, by = "gene")
# 
# saveRDS(genes, file = "genes.rds")
# genes <- readRDS(file = "genes.rds")

###################################################
### code chunk number : count matrix construct
###################################################
# library(GenomicAlignments)
# se <- summarizeOverlaps(features=genes, reads=bamfiles,mode = "Union",
#                        singleEnd=FALSE,ignore.strand=TRUE,fragments=TRUE)
# saveRDS(se, file = "se.rds")

###################################################
### code chunk number : banana counts data with meta data
###################################################
# se <- readRDS(file = "se.rds")
# colData(se) <- DataFrame(sampleTable)

###################################################
### code chunk number : DESeq data set
###################################################
# library(DESeq2)
# dds <- DESeqDataSet(se, design = ~ cell + day + condition)
# saveRDS(dds, file = "dds.rds")
dds <- readRDS(file = "dds.rds")

###################################################
### code chunk number : filter by counts
###################################################
# rs = rowSums(counts(dds))
# theta = 0.4
# use = (rs > quantile(rs, probs = theta))
# table(use)
# ddsFull = dds
# 
# quantile(rs,probs = 0.4)
# quantile(rs,probs = 0.2)
###################################################
### code chunk number : justify filter
###################################################
# ddsFull = DESeq(ddsFull)
# saveRDS(ddsFull, file = "ddsFull.rds")
# ddsFull = readRDS(file = "ddsFull.rds")
# pvalsFull = results(ddsFull)$pvalue
# plot(rank(rs)/length(rs),-log10(pvalsFull),pch=16,cex=0.45)
###################################################
### code chunk number : rlog for exploratory analysis
###################################################
# rld <- rlog(dds)
# saveRDS(rld, file = "rld.rds")
rld <- readRDS(file = "rld.rds")
# rld <- readRDS(file = "rldFilt.rds")
###################################################
### code chunk number : time points day6 and day8
###################################################
# idx_day6 = seq(1,23,2)
# idx_day8 = seq(2,24,2)

idx_day6 = seq(1,15,2)
idx_day8 = seq(2,16,2)
rld_day6 <- rld[, idx_day6]
rld_day8 <- rld[, idx_day8]

###################################################
### code chunk number : heatmap
###################################################
library("RColorBrewer")
library("gplots")
select = order(rowMeans(counts(dds)), decreasing=TRUE)[1:30]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)

###################################################
### code chunk number : sampleClust--sample distances
###################################################
mydist = function (x, method="euclidean")
{
        if (!is.na(pmatch(method, "euclidian"))) 
                method <- "euclidean"
        METHODS <- c("euclidean", "maximum", "manhattan", "canberra", 
                     "binary", "minkowski", "pearson", "spearman")
        method <- pmatch(method, METHODS)
        if (is.na(method)) 
                stop("invalid distance method")
        if (method == -1) 
                stop("ambiguous distance method")
        if (any(method == seq(1,6)))
        {
                dist(t(x), method = METHODS[method])
        }
        else
        {
                if (method == 7)
                        as.dist(1-cor(x, method = METHODS[method]))
                else
                        as.dist(1-cor(x, method = METHODS[method]))
        }           
}
rld_day6_assay = assay(rld_day6)
rld_day8_assay = assay(rld_day8)
rld_assay = assay(rld)

sampleDistsEuclidean_day6 <- as.matrix(mydist(rld_day6_assay))
sampleDistsEuclidean_day8 <- as.matrix(mydist(rld_day8_assay))
sampleDistsEuclidean <- as.matrix(mydist(rld_assay))

sampleDistsPearson_day6 <- as.matrix(mydist(rld_day6_assay,method = "pearson"))
sampleDistsPearson_day8 <- as.matrix(mydist(rld_day8_assay,method = "pearson"))
sampleDistsPearson <- as.matrix(mydist(rld_assay,method = "pearson"))

sampleDistsSpearman_day6 <- as.matrix(mydist(rld_day6_assay,method = "spearman"))
sampleDistsSpearman_day8 <- as.matrix(mydist(rld_day8_assay,method = "spearman"))
sampleDistsSpearman <- as.matrix(mydist(rld_assay,method = "spearman"))

# sampleDistsEuclidean_day6 <- as.matrix(dist(t(assay(rld_day6))))
# sampleDistsEuclidean_day8 <- as.matrix(dist(t(assay(rld_day8))))
# sampleDistsEuclidean <- as.matrix(dist(t(assay(rld))))
# 
# sampleDistsPearson_day6 <- 1-cor(assay(rld_day6), method = "pearson")
# sampleDistsPearson_day8 <- 1-cor(assay(rld_day8), method = "pearson")
# sampleDistsPearson <- 1-cor(assay(rld), method = "pearson")
# 
# sampleDistsSpearman_day6 <- 1-cor(assay(rld_day6), method =  "spearman")
# sampleDistsSpearman_day8 <- 1-cor(assay(rld_day8), method = "spearman")
# sampleDistsSpearman <- 1-cor(assay(rld), method = "spearman")

###################################################
### code chunk number : figHeatmapSamples--Euclidean
###################################################
# rownames(sampleDistsEuclidean_day6) = colnames(sampleDistsEuclidean_day6) = 
#         with(colData(rld_day6), paste(condition, cell, day, replicate, sep=":"))
# heatmap.2(sampleDistsEuclidean_day6, trace="none", col = rev(hmcol),
#           margin=c(13, 13),main = "Sample distances (Euclidean) day6")

# myheatmap = function (distMatrix, rld)
# {
#         var2str = function(v1)
#         {
#                 deparse(substitute(expr = v1))
#         }
#         rownames(distMatrix) = colnames(distMatrix) = 
#                 with(colData(rld), paste(cell, condition, day, replicate, sep=":"))
#         hc = hclust(as.dist(distMatrix), method = "average")
#         title = paste("Sample distances", var2str(distMatrix))
# #         title = paste("Sample distances", deparse(substitute(expr = distMatrix)))
#         print (title)
#         heatmap.2(distMatrix, Rowv= as.dendrogram(hc), symm=TRUE,
#                   trace="none", col = rev(hmcol), margin=c(13, 13),main = title)
# }
# myheatmap(sampleDistsEuclidean_day6, rld_day6)
rownames(sampleDistsEuclidean_day6) = colnames(sampleDistsEuclidean_day6) = 
  with(colData(rld_day6), paste(cell, condition, day, replicate, sep=":"))
hc = hclust(as.dist(sampleDistsEuclidean_day6), method = "average")
heatmap.2(sampleDistsEuclidean_day6, Rowv= as.dendrogram(hc), symm=TRUE,
          trace="none", col = rev(hmcol),
          margin=c(13, 13),main = "Sample distances (Euclidean) day6 NEW")

rownames(sampleDistsEuclidean_day8) = colnames(sampleDistsEuclidean_day8) = 
        with(colData(rld_day8), paste(cell, condition, day, replicate, sep=":"))
hc = hclust(as.dist(sampleDistsEuclidean_day8), method = "average")
heatmap.2(sampleDistsEuclidean_day8, Rowv= as.dendrogram(hc), symm=TRUE,
          trace="none", col = rev(hmcol),
          margin=c(13, 13),main = "Sample distances (Euclidean) day8 NEW")

rownames(sampleDistsEuclidean) = colnames(sampleDistsEuclidean) = 
        with(colData(rld), paste(cell, condition, day, replicate, sep=":"))
hc = hclust(as.dist(sampleDistsEuclidean), method = "average")
heatmap.2(sampleDistsEuclidean, Rowv= as.dendrogram(hc), symm=TRUE,
          trace="none", col = rev(hmcol),
          margin=c(13, 13),main = "Sample distances (Euclidean) full NEW")

# rownames(sampleDistsEuclidean_day6) = colnames(sampleDistsEuclidean_day6) = 
#         with(colData(rld_day6), paste(cell, condition, day, replicate, sep=":"))
# heatmap.2(sampleDistsEuclidean_day6, trace="none", col = rev(hmcol),
#           margin=c(13, 13),main = "Sample distances (Euclidean) day6")
# 
# rownames(sampleDistsEuclidean_day8) = colnames(sampleDistsEuclidean_day8) = 
#   with(colData(rld_day8), paste(cell, condition, day, replicate, sep=":"))
# heatmap.2(sampleDistsEuclidean_day8, trace="none", col = rev(hmcol),
#           margin=c(13, 13),main = "Sample distances (Euclidean) day8")
# 
# rownames(sampleDistsEuclidean) = colnames(sampleDistsEuclidean) = 
#   with(colData(rld), paste(cell, condition, day, replicate, sep=":"))
# heatmap.2(sampleDistsEuclidean, trace="none", col = rev(hmcol),
#           margin=c(13, 13),main = "Sample distances (Euclidean) full")

###################################################
### code chunk number : figHeatmapSamples--Pearson
###################################################
rownames(sampleDistsPearson_day6) = colnames(sampleDistsPearson_day6) = 
        with(colData(rld_day6), paste(cell, condition, day, replicate, sep=":"))
hc = hclust(as.dist(sampleDistsPearson_day6), method = "average")
heatmap.2(sampleDistsPearson_day6, Rowv= as.dendrogram(hc), symm=TRUE,
          trace="none", col = rev(hmcol),
          margin=c(13, 13),main = "Sample distances (Pearson) day6 NEW")

rownames(sampleDistsPearson_day8) = colnames(sampleDistsPearson_day8) = 
        with(colData(rld_day8), paste(cell, condition, day, replicate, sep=":"))
hc = hclust(as.dist(sampleDistsPearson_day8), method = "average")
heatmap.2(sampleDistsPearson_day8, Rowv= as.dendrogram(hc), symm=TRUE,
          trace="none", col = rev(hmcol),
          margin=c(13, 13),main = "Sample distances (Pearson) day8 NEW")

rownames(sampleDistsPearson) = colnames(sampleDistsPearson) = 
        with(colData(rld), paste(cell, condition, day, replicate, sep=":"))
hc = hclust(as.dist(sampleDistsPearson), method = "average")
heatmap.2(sampleDistsPearson, Rowv= as.dendrogram(hc), symm=TRUE,
          trace="none", col = rev(hmcol),
          margin=c(13, 13),main = "Sample distances (Pearson) full NEW")

# rownames(sampleDistsPearson_day6) = colnames(sampleDistsPearson_day6) = 
#   with(colData(rld_day6), paste(cell, condition, day, replicate, sep=":"))
# heatmap.2(sampleDistsPearson_day6, trace="none", col = rev(hmcol),
#           margin=c(13, 13),main = "Sample distances (Pearson) day6")
# 
# rownames(sampleDistsPearson_day8) = colnames(sampleDistsPearson_day8) = 
#   with(colData(rld_day8), paste(cell, condition, day, replicate, sep=":"))
# heatmap.2(sampleDistsPearson_day8, trace="none", col = rev(hmcol),
#           margin=c(13, 13),main = "Sample distances (Pearson) day8")
# 
# rownames(sampleDistsPearson) = colnames(sampleDistsPearson) = 
#   with(colData(rld), paste(cell, condition, day, replicate, sep=":"))
# heatmap.2(sampleDistsPearson, trace="none", col = rev(hmcol),
#           margin=c(13, 13),main = "Sample distances (Pearson) full")

###################################################
### code chunk number : figHeatmapSamples--Spearman
###################################################

rownames(sampleDistsSpearman_day6) = colnames(sampleDistsSpearman_day6) = 
        with(colData(rld_day6), paste(cell, condition, day, replicate, sep=":"))
hc = hclust(as.dist(sampleDistsSpearman_day6), method = "average")
heatmap.2(sampleDistsSpearman_day6, Rowv= as.dendrogram(hc), symm=TRUE,
          trace="none", col = rev(hmcol),
          margin=c(13, 13),main = "Sample distances (Spearman) day6 NEW")

rownames(sampleDistsSpearman_day8) = colnames(sampleDistsSpearman_day8) = 
        with(colData(rld_day8), paste(cell, condition, day, replicate, sep=":"))
hc = hclust(as.dist(sampleDistsSpearman_day8), method = "average")
heatmap.2(sampleDistsSpearman_day8, Rowv= as.dendrogram(hc), symm=TRUE,
          trace="none", col = rev(hmcol),
          margin=c(13, 13),main = "Sample distances (Spearman) day8 NEW")

rownames(sampleDistsSpearman) = colnames(sampleDistsSpearman) = 
        with(colData(rld), paste(cell, condition, day, replicate, sep=":"))
hc = hclust(as.dist(sampleDistsSpearman), method = "average")
heatmap.2(sampleDistsSpearman, Rowv= as.dendrogram(hc), symm=TRUE,
          trace="none", col = rev(hmcol),
          margin=c(13, 13),main = "Sample distances (Spearman) full NEW")
# rownames(sampleDistsSpearman_day6) = colnames(sampleDistsSpearman_day6) = 
#         with(colData(rld_day6), paste(cell, condition, day, replicate, sep=":"))
# heatmap.2(sampleDistsSpearman_day6, trace="none", col = rev(hmcol),
#           margin=c(13, 13),main = "Sample distances (Spearman) day6")
# 
# rownames(sampleDistsSpearman_day8) = colnames(sampleDistsSpearman_day8) = 
#         with(colData(rld_day8), paste(cell, condition, day, replicate, sep=":"))
# heatmap.2(sampleDistsSpearman_day8, trace="none", col = rev(hmcol),
#           margin=c(13, 13),main = "Sample distances (Spearman) day8")
# 
# rownames(sampleDistsSpearman) = colnames(sampleDistsSpearman) = 
#         with(colData(rld), paste(cell, condition, day, replicate, sep=":"))
# heatmap.2(sampleDistsSpearman, trace="none", col = rev(hmcol),
#           margin=c(13, 13),main = "Sample distances (Spearman) full")
###################################################
### code chunk number : figPCA
###################################################
library(ggplot2)
pca_day6 = plotPCA(rld_day6, intgroup=c("cell", "condition", "day", "replicate"), returnData = TRUE)
p = ggplot(data = pca_day6,
       aes_string(x = "PC1", y = "PC2", color = "cell", shape="condition") )
p = p + geom_point(size = 5)
p = p + xlab(paste0("PC1: ", round(attr(pca_day6, "percentVar")[1] * 100), "% variance")) + 
        ylab(paste0("PC2: ", round(attr(pca_day6, "percentVar")[2] * 100), "% variance"))
p = p + ggtitle(label = "PCA for day6")
p = p + theme(plot.title=element_text(face="bold",lineheight = 0.8))
p = p + geom_text(aes(label=paste(cell, condition, replicate, sep="")), vjust=0,hjust=0)
p

pca_day8 = plotPCA(rld_day8, intgroup=c("cell", "condition", "day", "replicate"), returnData = TRUE)
p = ggplot(data = pca_day8,
           aes_string(x = "PC1", y = "PC2", color = "cell", shape="condition") )
p = p + geom_point(size = 5)
p = p + xlab(paste0("PC1: ", round(attr(pca_day8, "percentVar")[1] * 100), "% variance")) + 
        ylab(paste0("PC2: ", round(attr(pca_day8, "percentVar")[2] * 100), "% variance"))
p = p + ggtitle(label = "PCA for day8")
p = p + theme(plot.title=element_text(face="bold",lineheight = 0.8))
p = p + geom_text(aes(label=paste(cell, condition, replicate, sep="")), vjust=0,hjust=0)
p

pca = plotPCA(rld, intgroup=c("cell", "condition", "day", "replicate"), returnData = TRUE)
p = ggplot(data = pca,
           aes_string(x = "PC1", y = "PC2", color = "cell", shape="condition") )
p = p + geom_point(size = 5)
p = p + xlab(paste0("PC1: ", round(attr(pca, "percentVar")[1] * 100), "% variance")) + 
        ylab(paste0("PC2: ", round(attr(pca, "percentVar")[2] * 100), "% variance"))
p = p + ggtitle(label = "PCA for all")
p = p + theme(plot.title=element_text(face="bold",lineheight = 0.8))
p = p + geom_text(aes(label=paste(cell, condition, replicate, sep="")), vjust=0,hjust=0)
p
###################################################
### code chunk number : MDS data
###################################################
library(scatterplot3d)
mds = data.frame(cmdscale(sampleDistsEuclidean,k = 3))
mds = cbind(mds, colData(rld))

mds_day6 = mds[idx_day6,]
mds_day8 = mds[idx_day8,]
# mds_day6 = mds[mds$day=="D6",]
# mds_day8 = mds[mds$day=="D8",]
s3d.data <- data.frame(dimension1=mds$X1, dimension2=mds$X2, dimension3=mds$X3)
s3d.data_day6 <- data.frame(dimension1=mds_day6$X1, dimension2=mds_day6$X2, dimension3=mds_day6$X3)
s3d.data_day8 <- data.frame(dimension1=mds_day8$X1, dimension2=mds_day8$X2, dimension3=mds_day8$X3)

###################################################
### code chunk number : figMDS
###################################################
library(colorspace)
mds_colors = rainbow_hcl(2)[c(mds_day6$cell)]
mds_pchs = c(mds_day6$condition)+15
s3d_day6 <- scatterplot3d(s3d.data_day6, type="h",color = mds_colors, lwd = 3, main="MDS plot for day6",pch="")
legend(s3d_day6$xyz.convert(0,0,max(mds_day6$X3)+10), col=unique(mds_colors),pch=16,legend = unique(mds_day6$cell),bty = "n")
text(s3d_day6$xyz.convert(s3d.data_day6), labels = mds_day6$replicate, offset=0.5,col = mds_colors,adj = c(0,0))
s3d_day6$points3d(s3d.data_day6, pch=mds_pchs, col = mds_colors, cex=1.5)
legend(s3d_day6$xyz.convert(-30,0,max(mds_day6$X3)+10), pch=mds_pchs[c(1,2)], legend=mds_day6$condition[c(1,2)],bty="n")

mds_colors = rainbow_hcl(2)[c(mds_day8$cell)]
mds_pchs = c(mds_day8$condition)+15
s3d_day8 <- scatterplot3d(s3d.data_day8, type="h",color = mds_colors, lwd = 3, main="MDS plot for day8",pch="")
legend(s3d_day8$xyz.convert(0,0,max(mds_day8$X3)+10), col=unique(mds_colors),pch=16,legend = unique(mds_day8$cell),bty = "n")
text(s3d_day8$xyz.convert(s3d.data_day8), labels = mds_day8$replicate, offset=0.5,col = mds_colors,adj = c(0,0))
s3d_day8$points3d(s3d.data_day8, pch=mds_pchs, col = mds_colors, cex=1.5)
legend(s3d_day8$xyz.convert(-30,0,max(mds_day8$X3)+10), pch=mds_pchs[c(1,2)], legend=mds_day8$condition[c(1,2)],bty="n")

mds_colors = rainbow_hcl(2)[c(mds$cell)]
mds_pchs = c(mds$condition)+15
s3d <- scatterplot3d(s3d.data, type="h",color = mds_colors, lwd = 3, main="MDS plot for all",pch="")
legend(s3d$xyz.convert(0,0,max(mds$X3)+10), col=unique(mds_colors),pch=16,legend = unique(mds$cell),bty = "n")
text(s3d$xyz.convert(s3d.data), labels = mds$replicate, offset=0.5,col = mds_colors,adj = c(0,0))
s3d$points3d(s3d.data, pch=mds_pchs, col = mds_colors, cex=1.5)
legend(s3d$xyz.convert(-30,0,max(mds$X3)+10), pch=mds_pchs[c(1,2)], legend=mds$condition[c(1,2)],bty="n")

###################################################
### code chunk number : dev.off
###################################################
dev.off()
# 
# mds <- data.frame(cmdscale(sampleDistMatrix, k=3))
# mds <- cbind(mds, colData(rld))
# mds <- mds[mds$day=="D6",]
# s3d.data <- data.frame(dimension1=mds$X1, dimension2=mds$X2, dimension3=mds$X3)
# library(colorspace)
# colors = rainbow_hcl(2)[c(mds$cell)]
# pchs = c(mds$condition)+15
# s3d <- scatterplot3d(s3d.data, type="h",color = colors, lwd = 3, main="MDS plot for day6",pch="")
# legend(s3d$xyz.convert(0,0,max(mds$X3)+10), col=unique(colors),pch=16,legend = unique(mds$cell),bty = "n")
# text(s3d$xyz.convert(s3d.data), labels = mds$replicate, offset=0.5,col = colors,adj = c(0,0))
# s3d$points3d(s3d.data, pch=pchs, col = colors, cex=1.5)
# legend(s3d$xyz.convert(-30,0,max(mds$X3)+10), pch=pchs[c(1,2)], legend=mds$condition[c(1,2)],bty="n")
# 
# mds <- data.frame(cmdscale(sampleDistMatrix, k=3))
# mds <- cbind(mds, colData(rld))
# mds <- mds[mds$day=="D8",]
# s3d.data <- data.frame(dimension1=mds$X1, dimension2=mds$X2, dimension3=mds$X3)
# library(colorspace)
# colors = rainbow_hcl(2)[c(mds$cell)]
# pchs = c(mds$condition)+15
# s3d <- scatterplot3d(s3d.data, type="h",color = colors, lwd = 3, main="MDS plot for day8",pch="")
# legend(s3d$xyz.convert(0,0,max(mds$X3)+35), col=unique(colors),pch=16,legend = unique(mds$cell),bty = "n")
# text(s3d$xyz.convert(s3d.data), labels = mds$replicate, offset=0.5,col = colors,adj = c(0,0))
# s3d$points3d(s3d.data, pch=pchs, col = colors, cex=1.5)
# legend(s3d$xyz.convert(-30,0,max(mds$X3)+35), pch=pchs[c(1,2)], legend=mds$condition[c(1,2)],bty="n")
# 
# mds <- data.frame(cmdscale(sampleDistMatrix, k=2))
# mds <- cbind(mds, colData(rld))
# # s2d.data <- data.frame(dimension1=mds$X1, dimension2=mds$X2)
# p <- qplot(X1, X2,  color=cell, shape=condition, data=mds, main="MDS plot-facets", size=I(5), facets=cell~day, label=replicate)
# p + geom_text(vjust=0,hjust=0)
