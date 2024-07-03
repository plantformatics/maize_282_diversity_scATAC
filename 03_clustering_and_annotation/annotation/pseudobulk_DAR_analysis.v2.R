## pseudo-bulk DEG analysis ##

# load libraries
library(SingleCellExperiment)
library(Matrix.utils)
library(edgeR)
library(parallel)

args <- commandArgs(T)
if(length(args) != 3){stop("Rscript pseudobulk_DAR_analysis.R <meta> <sparse.rds> <markers>")}

# functions
edgeRscDEG <- function(sce, ids=c("LouvainClusters", "library"), threads=1, use.ref.cells=F){
    
    # set up vars
    groups <- colData(sce)[, ids]
    groups[,1] <- as.factor(groups[,1])
    groups[,2] <- as.factor(groups[,2])
    clusters <- levels(groups[,1])
    threads <- ifelse(length(clusters) > threads, threads, length(clusters))
    
    # iterate
    outs <- mclapply(clusters, function(z){
        
        # set-up reference cells
        if(use.ref.cells){
            dff <- as.data.frame(groups)
            rownames(dff) <- rownames(colData(sce))
            not.cluster <- dff[dff[,1] != z,]
            cells.per.cluster <- nrow(dff[dff[,1] == z,])
            ran.cells <- sample(rownames(not.cluster), cells.per.cluster)
            #message(" - selected ", length(ran.cells), " control cells ...")
        }
        
        # verbose
        message(" - identifying DEG from cluster ", z)
        
        # rename groups
        df <- as.data.frame(groups)
        rownames(df) <- rownames(colData(sce))
        df[,1] <- as.character(df[,1])
        df$cluster_id <- as.factor(ifelse(df[,1] == z, 1, 0))
        df[,1] <- NULL
        if(use.ref.cells){
            df <- df[ifelse(df$cluster_id==1 | rownames(df) %in% ran.cells, T, F),]
            message(" - number of cluster cells = ", nrow(subset(df, df$cluster_id==1)))
            message(" - number of control cells = ", nrow(subset(df, df$cluster_id==0)))
        }
        #message(" - total nummber of cells in test = ", nrow(df))
        n.genes <- Matrix::colSums(counts(sce)[,rownames(df)] > 0)
        sample.data <- aggregate(n.genes~df$cluster_id+df[,ids[2]], FUN=mean)
        sample.data <- sample.data[order(sample.data[,1], decreasing=F),]
        
        # aggregate by cluster/library
        pb <- aggregate.Matrix(t(counts(sce)[,rownames(df)]), 
                               groupings = df, fun = "sum")
        pb <- t(pb)
        colnames(sample.data) <- c("cluster_id",ids[2], "n.genes")
        #print(head(pb))
        #print(head(sample.data))
        
        # create edgeR object
        group <- factor(sample.data$cluster_id, levels=c(0,1))
        batch <- factor(sample.data[,ids[2]])
        ave.genes <- scale(sample.data$n.genes)
        dge <- DGEList(counts = pb, 
                       norm.factors = rep(1, length(pb[1,])), 
                       group = group)

        # design experiment and estimate norm factors/dispersion
        design <- model.matrix(~ group)
        dge <- calcNormFactors(dge, method = "TMM", logratioTrim=0.1)
        dge <- estimateDisp(dge, design = design)
        
        # estimate Differential expression
        fit <- glmQLFit(dge, design)
        res <- glmQLFTest(fit, coef=ncol(fit$coefficients))
        res$table$FDR <- p.adjust(res$table[,4], method="fdr")
        res$table$cluster_id <- z
        res$table$geneID <- rownames(res$table)
        rownames(res$table) <- seq(1:nrow(res$table))
        message("   ~ returning results for cluster ",z)
        return(res$table)
        
    }, mc.cores=threads)
    outs <- do.call(rbind, outs)
    return(outs)
    
}

# load data
message(" - loading data ...")
counts <- readRDS(as.character(args[2]))
metadata <- read.table(as.character(args[1]))
shared <- intersect(rownames(metadata), colnames(counts))
metadata <- metadata[shared,]
counts <- counts[,shared]
counts <- counts[Matrix::rowSums(counts) > 0,]
sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata)
sce <- sce[rowSums(counts(sce) > 0) >= 0, ]
colData(sce)$batch <- colData(sce)$library
markers <- read.table(as.character(args[3]), header=T)
markers <- markers[!duplicated(markers$geneID),]
type <- markers$type
name <- markers$name
names(type) <- markers$geneID
names(name) <- markers$geneID

# filter low frequency, in-accessible genes
sce <- sce[rowMeans(counts(sce) > 0) > 0, ]

# iterate over all clusters
message(" - running DAR analysis by cluster with replicates ...")
deg <- edgeRscDEG(sce, threads=10, ids=c("leiden_refined", "library"), use.ref.cells=T)

# filter for significant genes (positive)
deg <- deg[order(deg$cluster_id, deg$PValue, decreasing=F),]
write.table(deg, file="DAR_pseudobulk.txt", quote=F, row.names=T, sep="\t", col.names=T)
sig <- subset(deg, deg$logFC > 0)
write.table(sig, file="DAR_pseudobulk.log2fc0.txt", quote=F, row.names=T, col.names=T, sep="\t")

# markers
sig.markers <- deg[as.character(deg$geneID) %in% as.character(markers$geneID),]
sig.markers$type <- type[sig.markers$geneID]
sig.markers$name <- name[sig.markers$geneID]
write.table(sig.markers, file="DAR_pseudobulk.marker_genes.txt", quote=F, row.names=T, col.names=T, sep="\t")
