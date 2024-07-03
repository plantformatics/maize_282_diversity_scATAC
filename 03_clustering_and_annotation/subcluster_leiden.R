## create shared nearest neighbor graph

# libraries
library(Seurat)
library(igraph)
library(uwot)
library(RColorBrewer)
library(mclust)
library(parallel)

# functions
doClusters <- function(svd, mdat, 
                       dims=ncol(svd), 
                       k=20, 
                       metric="cosine", 
                       res=2e-4, 
                       iters=100, 
                       modularity="CPM",
                       column="leiden_clusters"){
    
    faux <- Matrix(rep(0, length=100*nrow(svd)), ncol=nrow(svd), nrow=100, dimnames=list(seq(1:100), rownames(svd)))
    sro <- SeuratObject::CreateSeuratObject(faux, min.cells=0, min.features=0)
    sro[["svd"]] <- CreateDimReducObject(embeddings = as.matrix(svd), key = "PC_", assay = DefaultAssay(sro))
    
    # get snn
    sro <- FindNeighbors(sro, 
                         dims = 1:dims, 
                         reduction="svd", 
                         nn.eps=0, 
                         k.param=k, 
                         annoy.metric=metric,
                         n.trees=100,
                         prune.SNN=1/k,
                         l2.norm=F)
    
    # run with optimized resolution
    graph.obj <- graph_from_adjacency_matrix(as(sro[[paste0(DefaultAssay(sro),"_nn")]], "dgCMatrix"),
                                             mode="undirected",
                                             weighted=NULL)
    clusters <- cluster_leiden(graph.obj, 
                               objective_function=modularity,
                               resolution_parameter=res,
                               n_iterations=iters)
    clustered <- membership(clusters)
    mdat.f <- mdat
    mdat.f[,c(column)] <- clustered[rownames(mdat.f)]
    message(" - number of clusters = ", length(unique(mdat.f[,c(column)])))
    return(mdat.f)
}
refineClusters <- function(df, pca, method="svd", prop.z= -1){
    
    # define class
    df$leiden_clusters <- as.numeric(df$leiden_clusters)
    
    # set reduced dimensions
    if(method=="umap"){
        mat <- as.matrix(df[,c("umap1","umap2")])
    }else if(method=="svd"){
        mat <- as.matrix(pca)
    }
    
    # get cluster centroids
    clusts <- seq(1:max(df$leiden_clusters))
    ccs <- lapply(clusts, function(x){
        ids <- rownames(df[df$leiden_clusters==x,])
        if(length(ids) > 1){
            colMeans(mat[ids,])
        }else{
            as.numeric(mat[ids,])
        }
    })
    ccs <- do.call(rbind, ccs)
    rownames(ccs) <- paste0("cluster",clusts)
    
    # get counts per cluster
    num.cells.cluster <- prop.table(table(df$leiden_clusters))
    num.cells.cluster <- (num.cells.cluster - mean(num.cells.cluster))/sd(num.cells.cluster)
    merge.ids <- paste0("cluster",names(num.cells.cluster)[num.cells.cluster < prop.z])
    if(length(num.cells.cluster[num.cells.cluster < prop.z]) < 1){
        return(df)
    }else{
    
        # cluster-cluster distance
        cc.dist <- as.matrix(dist(ccs, diag=T, upper=T))
        
        # append clusters with low cell counts to nearest cluster
        matched <- cc.dist[merge.ids,!colnames(cc.dist) %in% merge.ids]
        matched <- apply(matched, 1, function(x){names(x)[which.min(x)]})
        for(i in 1:length(matched)){
            id <- names(matched)[i]
            id <- as.numeric(gsub("cluster","",id))
            new <- matched[i]
            new <- as.numeric(gsub("cluster","",new))
            df$leiden_clusters <- ifelse(df$leiden_clusters==id, new, df$leiden_clusters)
        }
        
        # update cluster IDs
        new.cids <- table(df$leiden_clusters)
        new.cids <- new.cids[order(new.cids, decreasing=T)]
        names(new.cids) <- paste0("cluster", names(new.cids))
        df$leiden_clusters2 <- factor(paste0("cluster",df$leiden_clusters), levels=names(new.cids))
        df$leiden_clusters2 <- as.numeric(df$leiden_clusters2)
        df$leiden_clusters <- df$leiden_clusters2
        df$leiden_clusters2 <- NULL
        
        # return
        return(df)
    }
}
pruneClusters <- function(df, threshold=5, threads=10){
    
    # filter outliers?
    .filterSingle  <- function(pro, k=30, threshold=3){
        
        # set column names
        vars <- colnames(pro)
        
        # ensure that k is less than number of samples
        if(k > nrow(pro)){
            k <- nrow(pro)-1
        }
        
        # get nearest neighbors
        topk <- FNN::get.knn(pro[,vars], k=k)
        cell.dists <- as.matrix(topk$nn.dist)
        rownames(cell.dists) <- rownames(pro)
        colnames(cell.dists) <- paste0("k",seq(1:ncol(cell.dists)))
        zscore <- (cell.dists - mean(cell.dists))/sd(cell.dists)
        names(zscore) <- rownames(pro)
        ave <- rowMeans(zscore)
        
        # thresholds
        p.zscore <- ave[order(ave, decreasing=T)]
        num.pass <- length(ave[ave < threshold])
        
        # filter
        prop.good <- ave[ave < threshold]
        ids <- names(prop.good)
        
        # return object
        return(ids)
    }
    
    # iterate over clusters
    clusters <- seq(1:max(df$leiden_clusters))
    ids <- mclapply(clusters, function(x){
        message(" - pruning cells in cluster ",x, "/",max(df$leiden_clusters))
        df.sub <- df[df$leiden_clusters==x,]
        .filterSingle(df[,c("umap1","umap2")], k=30, threshold=threshold)
    }, mc.cores=threads)
    ids <- unlist(ids)
    df$prune_filter <- ifelse(rownames(df) %in% ids, 1, 0)
    return(df)
}
findRes <- function(svd, 
                    meta,
                    metric="cosine",
                    k=30,
                    res=seq(from=1e-4,to=5e-3, length.out=100), 
                    perms=5){
    
    # get NN graph
    faux <- Matrix(rep(0, length=100*nrow(svd)), ncol=nrow(svd), nrow=100, dimnames=list(seq(1:100), rownames(svd)))
    sro <- SeuratObject::CreateSeuratObject(faux, min.cells=0, min.features=0)
    sro[["svd"]] <- CreateDimReducObject(embeddings = as.matrix(svd), key = "PC_", assay = DefaultAssay(sro))
    
    # get snn
    sro <- FindNeighbors(sro, 
                         dims = 1:ncol(svd), 
                         reduction="svd", 
                         nn.eps=0, 
                         k.param=k, 
                         annoy.metric=metric,
                         n.trees=100,
                         prune.SNN=1/k,
                         l2.norm=F)
    
    # make graph object
    graph.obj <- graph_from_adjacency_matrix(as(sro[[paste0(DefaultAssay(sro),"_nn")]], "dgCMatrix"),
                                             mode="undirected",
                                             weighted=NULL)
    
    # prep perturbed graph
    perturb.graph <- as(sro[[paste0(DefaultAssay(sro),"_nn")]], "dgCMatrix")
    p.g.df.in <- as.data.frame(summary(perturb.graph))
    
    # iterate over different resolutions
    outs <- mclapply(res, function(z){
        
        # message
        message(' - checking resolution = ', z)
        
        # perturb 2% of edges over x iterations
        p.graphs <- lapply(seq(1:perms), function(xx){
            set.seed((xx+1111))
            p.g.df <- p.g.df.in[p.g.df.in$i != p.g.df.in$j,]
            num.rows <- ceiling(nrow(p.g.df)*0.02)
            p.g.df$x[sample(nrow(p.g.df), num.rows)] <- 0
            p.g.df.0 <- subset(p.g.df, p.g.df$x == 0)
            p.g.df.1 <- data.frame(i=p.g.df.0$j, j=p.g.df.0$i, x=p.g.df.0$x)
            p.g.df.new <- rbind(p.g.df.0, p.g.df.1, p.g.df)
            p.g.df.new <- p.g.df.new[order(p.g.df.new$x, decreasing=F),]
            p.g.df.new <- p.g.df.new[!duplicated(p.g.df.new[,1:2]),]
            perturb.graph2 <- sparseMatrix(i=p.g.df.new$i,
                                           j=p.g.df.new$j,
                                           x=p.g.df.new$x,
                                           dimnames=list(rownames(sro[["svd"]]),
                                                         rownames(sro[["svd"]])))
            graph_from_adjacency_matrix(as(perturb.graph2, "dgCMatrix"),
                                        mode="undirected",
                                        weighted=NULL)
        })
        
        # cluster
        graph.obj.1 <- cluster_leiden(graph.obj, 
                                      objective_function="CPM",
                                      resolution_parameter=z,
                                      n_iterations=2)
        clustered <- membership(graph.obj.1)
        message(" - # clusters = ", length(table(clustered)))
        num.clust <- length(table(clustered))
        
        # cluster perturbed
        p.clusts <- lapply(p.graphs, function(xx){
            perturb.1 <- cluster_leiden(xx,
                                        objective_function="CPM",
                                        resolution_parameter=z,
                                        n_iterations=2)
            as.numeric(as.factor(membership(perturb.1)))
        })
        
        # calc average ARI
        message("   calculating ARI ...")
        ari <- vector(length=length(p.clusts))
        for(i in 1:length(p.clusts)){
            ari[i] <- adjustedRandIndex(clustered, p.clusts[[i]])
        }
        ave.ari <- mean(ari)
        data.frame(res=z, ari=ave.ari, total.clusters=num.clust)
        
    }, mc.cores=30)
    outs <- do.call(rbind, outs)
    return(outs)
    
}
clusterConsensus <- function(sro, 
                             res=seq(from=1e-4,to=5e-3, length.out=100), 
                             rounds=100, 
                             threads=30){
    
    # make graph object
    graph.obj <- graph_from_adjacency_matrix(as(sro[[paste0(DefaultAssay(sro),"_nn")]], "dgCMatrix"),
                                             mode="undirected",
                                             weighted=NULL)
    # iterate over different resolutions
    outs <- mclapply(res, function(z){
        
        # message
        message(' - checking resolution = ', z)
        
        # create connectivity matrix
        con <- Matrix(rep(0, length=nrow(sro[["svd"]])**2), 
                      ncol=nrow(sro[["svd"]]), 
                      nrow=nrow(sro[["svd"]]), 
                      dimnames=list(rownames(sro[["svd"]]), rownames(sro[["svd"]])))
        
        # iterate over various seeds
        for (x in seq(1:rounds)){
            
            # verbose
            message("   round ", x)
            # set seed
            set.seed(x)
            
            # cluster
            graph.obj.1 <- cluster_leiden(graph.obj, 
                                          objective_function="CPM",
                                          resolution_parameter=z,
                                          n_iterations=10)
            clustered <- membership(graph.obj.1)
            for(i in 1:length(clustered)){
                for(j in 1:length(clustered)){
                    if(clustered[i]==clustered[j]){
                        con[names(clustered)[i],names(clustered)[j]] <- con[names(clustered)[i],names(clustered)[j]] + 1
                    }
                }
            }
        }
        return(con/100)
        
    }, mc.cores=threads)
    return(outs)
}
callClusterACRs <- function(bed, metadata, cluster.var="leiden_clusters"){
    
    # ACR calling function
    .callACRs <- function(obj, genomesize=1.6e9, shift= -50, extsize=100,
                         output="bulk_peaks", tempdir="./macs2_temp", verbose=T, fdr=0.05){
        
        # verbose
        if(verbose){message(" - running MACS2 on bulk BED file ...")}
        bed <- obj$bedpath
        
        # create temp dif
        mac2temp <- tempdir
        if(file.exists(mac2temp)){
            unlink(mac2temp)
            dir.create(mac2temp)
        }else{
            dir.create(mac2temp)
        }
        
        # build command
        cmdline <- paste0("macs2 callpeak -t ", bed, " -f BED -g ", genomesize, " --keep-dup all -n ", output,
                          " --nomodel --shift ",shift, " --extsize ", extsize, " --outdir ",mac2temp, " --qvalue ", fdr)
        
        # run macs2
        suppressMessages(system(cmdline))
        
        # load peaks
        peaks <- read.table(paste0(mac2temp,"/",output,"_peaks.narrowPeak"))
        
        # delete temp dir
        unlink(mac2temp)
        
        # return acrs
        return(peaks)
        
    }
    
    # load bed file
    obj <- list()
    if(grepl(".gz$", bed)){
        obj$bed <- read.table(gzfile(as.character(bed)))
    }else{
        obj$bed <- read.table(as.character(bed))
    }
    
    # call peaks per cluster
    peaks <- lapply(unique(metadata[,c(cluster.var)]), function(z){
        ids <- rownames(metadata[metadata[,c(cluster.var)]==z,])
        cl.bed <- obj$bed[as.character(obj$bed$V4) %in% ids,]
        
    })
    
    # find overlaps in peaks
    message(" - building sparse ACR matrix ...")
    acrs <- GRanges(seqnames=obj$acr$V1,
                    ranges=IRanges(start=obj$acr$V2,
                                   end=obj$acr$V3,
                                   names=paste(obj$acr$V1,obj$acr$V2,obj$acr$V3, sep="_")))
    peak_hits <- suppressWarnings(findOverlaps(Tn5_Grange, acrs, 
                                               minoverlap=1, 
                                               type=c("within"), 
                                               select="all", 
                                               ignore.strand=T))
    
    peak_int <- paste(names(acrs)[peak_hits@to],
                      names(Tn5_Grange)[peak_hits@from],sep="/_Com_/")
    peak_int <- table(peak_int)
    peak_int <- data.frame(acrID = as.factor(as.character(lapply(strsplit(as.character(rownames(peak_int)),
                                                                          split="/_Com_/"), "[", 1))),
                           cellID = as.factor(as.character(lapply(strsplit(as.character(rownames(peak_int)),
                                                                           split="/_Com_/"), "[", 2))),
                           activity = as.numeric(as.character(peak_int)))
    obj$pa <- sparseMatrix(i=as.numeric(peak_int$acrID),
                           j=as.numeric(peak_int$cellID),
                           x=peak_int$activity,
                           dimnames=list(levels(peak_int$acrID),
                                         levels(peak_int$cellID)))
    
}
clstcolors <- function(cols){
    
    if(cols=="grove"){
        return(c("#017351", "#fbbf45", "#ed0345", "#3B9AB2", "#710162", "#03c383", 
                 "#a12a5e", "#aad962", "#01545a", "#ef6a32", "#1a1334"))
    }
    if(cols=="bear"){
        return(c("#faa818", "#41a30d", "#fbdf72", "#367d7d", "#d33502", "#6ebcbc", 
                 "#37526d", "#916848", "#f5b390", "#342739", "#bed678", "#a6d9ee", 
                 "#0d74b6", "#60824f", "#725ca5", "#e0598b"))
    }
    
}
shadowtext <- function(x, y=NULL, 
                       labels, 
                       col='white', bg='black', 
                       theta= seq(0, 2*pi, length.out=50), 
                       r=0.1, ... ) {
    
    xy <- xy.coords(x,y)
    xo <- r*strwidth('A')
    yo <- r*strheight('A')
    
    # draw background text with small shift in x and y in background colour
    for (i in theta) {
        text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
    }
    # draw actual text in exact xy position in foreground colour
    text(xy$x, xy$y, labels, col=col, ... )
}

# load data
pcs <- read.table("maize_282.v8.3.ALL_CELLs.reduced_dimensions.txt")
meta <- read.table("maize_282.v8.3.ALL_CELLs.metadata.leiden.pruned.txt")
pcs.og <- pcs
pcs <- pcs.og

# variables
k <- 30
metric <- "cosine"
dims <- 20

# subset pcs
pcs <- pcs[,1:dims]
pcs <- t(apply(pcs, 1, function(x){(x-mean(x))/sd(x)}))

# match ids
shared <- intersect(rownames(pcs), rownames(meta))
pcs <- pcs[shared,]
meta <- meta[shared,]

# umap
ucor <- umap(pcs, metric=metric, n_neighbors=k, a=1.95, b=0.7, verbose=T)

# update columns
meta$umap1_og.1 <- meta$umap1
meta$umap2_og.1 <- meta$umap2
meta$umap1 <- ucor[rownames(meta),1]
meta$umap2 <- ucor[rownames(meta),2]

# get optimized resolution
#res <- findRes(pcs, meta, metric=metric, k=k, res=seq(from=0,to=5e-3,length.out=100))

# set up data
meta.f <- doClusters(pcs, meta, dims=dims, k=k, metric=metric, res=2.5, modularity="modularity")
meta.f <- refineClusters(meta.f, pcs, method="svd", prop.z = -1)
meta.ff <- pruneClusters(meta.f, threshold=3)

# save output
write.table(meta.ff, file="maize_282.v8.3.ALL_CELLs.metadata.leiden.pruned.txt")

# remove "pruned" cells
meta.f <- subset(meta.ff, meta.ff$prune_filter==1)

# plot prune
pdf("maize_282.v8.3.ALL_CELLs.LEIDEN.pruned.pdf", width=10, height=10)
plot(meta.ff$umap1, meta.ff$umap2, pch=16, cex=0.2, col=ifelse(meta.ff$prune_filter==1, "grey75", "black"))
dev.off()

# plot attributes
cl.coord <- lapply(unique(meta.f$leiden_clusters), function(x){
    df <- subset(meta.f, meta.f$leiden_clusters==x)
    x1 <- median(df$umap1, na.rm=T)
    y1 <- median(df$umap2, na.rm=T)
    out <- data.frame(x=x1,y=y1,cluster=x)
    return(out)
})
cl.coord <- do.call(rbind, cl.coord)

# plot
pdf("maize_282.v8.3.ALL_CELLs.leiden_clustering.pdf", width=10, height=10)
cols <- colorRampPalette(clstcolors("grove"))(length(unique(meta.f$leiden_clusters)))
x.max <- max(meta.f$umap1)
x.min <- min(meta.f$umap1)
plot(meta.f$umap1, meta.f$umap2, col=cols[meta.f$leiden_clusters], pch=16, cex=0.2,
     xlim=c(x.min, (x.max*1.5)))
shadowtext(x=cl.coord$x, y=cl.coord$y, labels=cl.coord$cluster, col="black", bg="white", cex=0.8)
legend("right", legend=seq(1:max(meta.f$leiden_clusters)), fill=cols, border=NA, col=cols)
dev.off()

# subcluster 
reclust <- c(1, 2, 5, 12, 
             13, 15, 27, 38)
resolutions <- c(0.3, 0.2, 0.2, 0.3,
                 0.35, 0.4, 0.25, 0.4)
j <- 0
pdf("subclusters_leiden_scATAC.pdf", width=16, height=8)
layout(matrix(c(1:length(reclust)), nrow=2, byrow=F))
new <- lapply(reclust, function(i){
    
    # select cluster
    j <<- j + 1
    meta.f.sub <- subset(meta.f, meta.f$leiden_clusters==i)
    pcs.f.sub <- pcs[rownames(meta.f.sub),]
    meta.f.sub <- doClusters(pcs.f.sub, meta.f.sub, res=resolutions[j], modularity="modularity", column="subcluster")
    plot(meta.f.sub$umap1, meta.f.sub$umap2, pch=16, cex=0.2, col=meta.f.sub$subcluster, main=i)
    return(meta.f.sub)
})
dev.off()

new <- do.call(rbind, new)
new$leiden_refined <- paste0(new$leiden_cluster, ".", new$subcluster)
clustID <- new$leiden_refined
names(clustID) <- rownames(new)
meta.f$leiden_refined <- clustID[rownames(meta.f)]
meta.f$leiden_refined <- ifelse(is.na(meta.f$leiden_refined), paste0(meta.f$leiden_clusters,".",0), meta.f$leiden_refined)

# plot attributes
meta.f <- subset(meta.f, meta.f$prune_filter==1)
cl.coord <- lapply(unique(meta.f$leiden_refined), function(x){
    df <- subset(meta.f, meta.f$leiden_refined==x)
    x1 <- median(df$umap1, na.rm=T)
    y1 <- median(df$umap2, na.rm=T)
    out <- data.frame(x=x1,y=y1,cluster=x)
    return(out)
})
cl.coord <- do.call(rbind, cl.coord)

# plot subclustered
pdf("maize_282.v8.3.ALL_CELLs.leiden_SUBclustering.pdf", width=10, height=10)
meta.f$leiden_refined <- as.factor(meta.f$leiden_refined)
cols <- colorRampPalette(clstcolors("grove"))(length(unique(meta.f$leiden_refined)))
names(cols) <- sort(unique(meta.f$leiden_refined))
x.max <- max(meta.f$umap1)
x.min <- min(meta.f$umap1)
plot(meta.f$umap1, meta.f$umap2, col=cols[meta.f$leiden_refined], pch=16, cex=0.2,
     xlim=c(x.min, (x.max*1.5)))
shadowtext(x=cl.coord$x, y=cl.coord$y, labels=cl.coord$cluster, col="black", bg="white", cex=0.7)
legend("right", legend=names(cols), fill=cols, border=NA, col=cols, cex=0.8)
dev.off()
