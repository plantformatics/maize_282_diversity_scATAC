## plot marker heatmap ##

# load libraries
library(pheatmap)
library(RColorBrewer)
library(edgeR)
library(Matrix)

# load arguments
args <- commandArgs(T)
if(length(args)!=3){stop("Rscript PlotClusterZscore.R <metadata> <gene.sparse.rds> <markers>")}

# load data
a <- read.table(as.character(args[1]))
markers <- read.table(as.character(args[3]), header=T)
b <- readRDS(as.character(args[2]))

# aggregate clusters
clusts <- sort(unique(a$leiden_refined))
out <- lapply(clusts, function(x){
    df <- subset(a, a$leiden_refined==x)
    Matrix::rowSums(b[,rownames(df)])
})
out <- do.call(cbind, out)
colnames(out) <- paste0("cluster_", clusts)

# normalize
out <- cpm(out, log=F)
gene.aves <- rowMeans(out)
out <- out[gene.aves > 0,]
vars <- apply(out, 1, var)
out <- out[vars > 0,]

# estimate log2FC
ids <- colnames(out)
fc <- lapply(ids, function(x){
    log2((out[,x])/rowMeans(out[,!colnames(out) %in% x]))
})
fc <- do.call(cbind, fc)
colnames(fc) <- ids
fc[is.infinite(fc) & fc < 0] <- -1*max(fc[is.finite(fc)])
fc[is.infinite(fc) & fc > 0] <- max(fc[is.finite(fc)])
fc[is.na(fc)] <- 0

# get zscores
zscore <- as.matrix(t(scale(t(out))))

# save tables
write.table(out, file="maize_282.v8.3_leiden_clusters.CPM.txt", quote=F, row.names=T, col.names=T, sep="\t")
write.table(fc, file="maize_282.v8.3_leiden_clusters.log2FC.txt", quote=F, row.names=T, col.names=T, sep="\t")
write.table(zscore, file="maize_282.v8.3_leiden_clusters.zscore.txt", quote=F, row.names=T, col.names=T, sep="\t")

# cluster columns
fc.clust <- hclust(dist(t(fc)))$order
z.clust <- hclust(dist(t(zscore)))$order
fc <- fc[,fc.clust]
zscore <- zscore[,z.clust]

# subset to markers
zscore.m <- zscore[rownames(zscore) %in% as.character(markers$geneID),]
fc.m <- fc[rownames(fc) %in% as.character(markers$geneID),]

# cluster rows
f.row <- apply(fc.m, 1, which.max)
z.row <- apply(zscore.m, 1, which.max)
fc.m <- fc.m[order(f.row, decreasing=F),]
zscore.m <- zscore.m[order(z.row, decreasing=F),]

# geneIDs
gname <- markers$name
names(gname) <- markers$geneID
rownames(fc.m) <- gname[rownames(fc.m)]
rownames(zscore.m) <- gname[rownames(zscore.m)]

# plot
pdf("Log2FC_marker_genes.heatmap.pdf", width=6, height=12)
pheatmap(fc.m, col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         cluster_col=F, cluster_row=F)
dev.off()

pdf("Zscore_marker_genes.heatmap.pdf", width=6, height=12)
val <- max(abs(zscore.m))
pheatmap(zscore.m, col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         cluster_col=F, cluster_row=F, breaks=seq(from= -val, to=val, length.out=101),
         fontsize_row=8)
dev.off()




