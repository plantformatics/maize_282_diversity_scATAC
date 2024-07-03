# run Socrates on merged socrates object #

# libraries
library(Socrates)
library(harmony)
library(symphony)
library(igraph)
library(Matrix)


# functions --------------------------------------------------------------

# remove high conf cells form raw data set
cleanObj <- function(raw, 
                     conf){
  
  # identify shared cells
  shared <- intersect(rownames(raw$meta), rownames(conf$Clusters))
  raw$counts <- raw$counts[,!colnames(raw$counts) %in% shared]
  shared.sites <- intersect(rownames(raw$counts), rownames(conf$residuals))
  raw$counts <- raw$counts[shared.sites,]
  raw$counts <- raw$counts[,Matrix::colSums(raw$counts)>0]
  raw$counts <- raw$counts[Matrix::rowSums(raw$counts)>0,]
  raw$meta <- raw$meta[colnames(raw$counts),]
  
  # return
  return(raw)
  
}

# extract feature information for symphony integration
extractFeatureMetrics <- function(obj, 
                                  num.var=0){
  
  # internal functions
  RowVar <- function(x) {
    spm <- t(x)
    stopifnot(methods::is(spm, "dgCMatrix"))
    ans <- sapply(base::seq.int(spm@Dim[2]), function(j) {
      if (spm@p[j + 1] == spm@p[j]) {
        return(0)
      }
      mean <- base::sum(spm@x[(spm@p[j] + 1):spm@p[j +
                                                     1]])/spm@Dim[1]
      sum((spm@x[(spm@p[j] + 1):spm@p[j + 1]] - mean)^2) +
        mean^2 * (spm@Dim[1] - (spm@p[j + 1] - spm@p[j]))
    })/(spm@Dim[1] - 1)
    names(ans) <- spm@Dimnames[[2]]
    ans
  }
  
  # get sites used for clustering
  residuals_slotName <- "residuals"
  row.var <- RowVar(obj[[residuals_slotName]])
  row.means <- Matrix::rowMeans(obj[[residuals_slotName]])
  adj.row.var <- loess(row.var ~ row.means)$residuals
  names(adj.row.var) <- names(row.var)
  adj.row.var <- adj.row.var[order(adj.row.var, decreasing = T)]
  if(num.var < 100){
    topSites <- names(adj.row.var[adj.row.var > num.var])
  }else{
    topSites <- names(adj.row.var)[1:num.var]
  }
  features <- data.frame(symbol=topSites, mean=row.means[topSites], stddev=sqrt(row.var[topSites]))
  loadings <- obj$PCA_model$v
  colnames(loadings) <- paste0("PC_",seq(1:ncol(loadings)))
  loadings <- loadings[,obj$PCA_model$keep_pcs]
  rownames(loadings) <- obj$hv_sites
  pvars <- obj$PCA_model$d[obj$PCA_model$keep_pcs]
  
  # return
  return(list(features=features, loadings=loadings, pvars=pvars))
  
}

# new TFIDF methods
tfidf <- function(obj,
                  frequencies=T,
                  log_scale_tf=T,
                  scale_factor=10000,
                  doL2=F,
                  slotName="residuals"){
  
  # set bmat
  bmat <- obj$counts
  
  # hidden functions
  .safe_tfidf       <- function(tf, idf,  block_size=2000e6){
    result = tryCatch({
      result = tf * idf
      result
    }, error = function(e) {
      options(DelayedArray.block.size=block_size)
      DelayedArray:::set_verbose_block_processing(TRUE)
      
      tf = DelayedArray(tf)
      idf = as.matrix(idf)
      
      result = tf * idf
      result
    })
    return(result)
  }
  
  # Use either raw counts or divide by total counts in each cell
  if (frequencies) {
    # "term frequency" method
    tf = t(t(bmat) / Matrix::colSums(bmat))
  } else {
    # "raw count" method
    tf = bmat
  }
  
  # Either TF method can optionally be log scaled
  if (log_scale_tf) {
    if (frequencies) {
      tf@x = log1p(tf@x * scale_factor)
    } else {
      tf@x = log1p(tf@x * 1)
    }
  }
  
  # IDF w/ "inverse document frequency smooth" method
  idf = log(1 + ncol(bmat) / Matrix::rowSums(bmat))
  
  # TF-IDF
  tf_idf_counts = .safe_tfidf(tf, idf)
  
  # do L2?
  if(doL2){
    l2norm <- function(x){x/sqrt(sum(x^2))}
    colNorm <- sqrt(Matrix::colSums(tf_idf_counts^2))
    tf_idf_counts <- tf_idf_counts %*% Diagonal(x=1/colNorm)
  }
  
  rownames(tf_idf_counts) = rownames(bmat)
  colnames(tf_idf_counts) = colnames(bmat)
  obj[[slotName]] <- Matrix(tf_idf_counts, sparse=T)
  obj$norm_method <- "tfidf"
  
  # return
  return(obj)
}

# use ref cells IDF to project new data
projectTFIDF <- function(obj, 
                         old.obj, 
                         doL2=T, 
                         slotName="residuals", 
                         scale_factor=10000){
  
  # set bmat
  bmat <- obj$counts
  old.obj$counts <- old.obj$counts[intersect(rownames(old.obj$residuals), rownames(bmat)),]
  bmat <- bmat[rownames(old.obj$counts),]
  
  # hidden functions
  .safe_tfidf <- function(tf, idf,  block_size=2000e6){
    result = tryCatch({
      result = tf * idf
      result
    }, error = function(e) {
      options(DelayedArray.block.size=block_size)
      DelayedArray:::set_verbose_block_processing(TRUE)
      
      tf = DelayedArray(tf)
      idf = as.matrix(idf)
      
      result = tf * idf
      result
    })
    return(result)
  }
  
  # "term frequency" method
  tf = t(t(bmat) / Matrix::colSums(bmat))
  
  # TF log scaled
  tf@x = log1p(tf@x * scale_factor)
  
  # IDF w/ "inverse document frequency smooth" method
  idf = log(1 + ncol(old.obj$counts) / Matrix::rowSums(old.obj$counts))
  
  # TF-IDF
  tf_idf_counts = .safe_tfidf(tf, idf)
  
  # do L2?
  if(doL2){
    colNorm <- sqrt(Matrix::colSums(tf_idf_counts^2))
    tf_idf_counts <- tf_idf_counts %*% Diagonal(x=1/colNorm)
  }
  
  rownames(tf_idf_counts) = rownames(bmat)
  colnames(tf_idf_counts) = colnames(bmat)
  obj[[slotName]] <- Matrix(tf_idf_counts, sparse=T)
  obj$norm_method <- "tfidf"
  
  # return
  return(obj)
  
}

# new map query
mapQuery2 <- function(exp_query, 
                      metadata_query, 
                      ref_obj,
                      svd_diag = NULL,
                      vars = NULL, 
                      verbose = TRUE,
                      do_umap = TRUE, 
                      sigma = 0.1,
                      scaleVar = T,
                      doSTD = F,
                      doL2 = T){
  
  # hidden functions
  cosine_normalize_cpp <- function(V, dim) {
    .Call('_symphony_cosine_normalize_cpp', PACKAGE = 'symphony', V, dim)
  }
  soft_cluster <- function(Y, Z, sigma) {
    .Call('_symphony_soft_cluster', PACKAGE = 'symphony', Y, Z, sigma)
  }
  moe_correct_ref <- function(Zq, Xq, Rq, Nr, RrZtr) {
    .Call('_symphony_moe_correct_ref', PACKAGE = 'symphony', Zq, Xq, Rq, Nr, RrZtr)
  }
  
  # verbose
  if (verbose){message("Scaling and synchronizing query gene expression")}
  
  # find shared features
  idx_shared_genes <- which(ref_obj$vargenes$symbol %in% rownames(exp_query))
  shared_genes <- ref_obj$vargenes$symbol[idx_shared_genes]
  
  # verbose
  if (verbose){
    message("Found ", length(shared_genes), " out of ", length(ref_obj$vargenes$symbol),
            " reference variable genes in query dataset")
  }
  
  # align matrices
  exp_query_scaled_sync <- exp_query[shared_genes,]
  rownames(exp_query_scaled_sync) <- ref_obj$vargenes$symbol
  colnames(exp_query_scaled_sync) <- colnames(exp_query)
  
  # verbose
  if (verbose){message("Project query cells using reference gene loadings")}
  
  # scale TF-IDF by feature loadings
  Z_pca_query <- as.matrix(t(ref_obj$loadings[shared_genes,]) %*% exp_query_scaled_sync)
  
  # scale PCs
  if(scaleVar & !is.null(svd_diag)){
    if (verbose){message("Scaling PCs by variance")}
    Z_pca_query <- t(t(Z_pca_query) %*% diag(svd_diag))
  }
  if(doSTD){
    Z_pca_query <- apply(Z_pca_query, 2, function(xx){(xx-mean(xx, na.rm=T))/sd(xx, na.rm=T)})
  }
  if(doL2){
    Z_pca_query <- apply(Z_pca_query, 2, function(xx){(xx/sqrt(sum(xx^2, na.rm=T)))})
  }
  
  # verbose
  if (verbose){message("Clustering query cells to reference centroids")}
  
  # normalize cosine distances
  Z_pca_query_cos <- cosine_normalize_cpp(Z_pca_query, 2)
  R_query <- soft_cluster(ref_obj$centroids, Z_pca_query_cos, sigma)
  
  # verbose
  if (verbose){message("Correcting query batch effects")}
  if (!is.null(vars)) {
    design <- droplevels(metadata_query)[, vars] %>% as.data.frame()
    onehot <- design %>% purrr::map(function(.x) {
      if (length(unique(.x)) == 1) {
        rep(1, length(.x))
      }
      else {
        stats::model.matrix(~0 + .x)
      }
    }) %>% purrr::reduce(cbind)
    Xq = cbind(1, intercept = onehot) %>% t()
  }
  else {
    Xq = Matrix(rbind(rep(1, ncol(Z_pca_query)), rep(1, ncol(Z_pca_query))),
                sparse = TRUE)
  }
  Zq_corr = moe_correct_ref(as.matrix(Z_pca_query), as.matrix(Xq),
                            as.matrix(R_query), as.matrix(ref_obj$cache[[1]]), as.matrix(ref_obj$cache[[2]]))
  colnames(Z_pca_query) = row.names(metadata_query)
  rownames(Z_pca_query) = paste0("PC_", seq_len(nrow(Zq_corr)))
  colnames(Zq_corr) = row.names(metadata_query)
  rownames(Zq_corr) = paste0("harmony_", seq_len(nrow(Zq_corr)))
  umap_query = NULL
  if (do_umap & !is.null(ref_obj$save_uwot_path)) {
    if (verbose)
      message("UMAP")
    ref_umap_model = uwot::load_uwot(ref_obj$save_uwot_path,
                                     verbose = FALSE)
    umap_query = uwot::umap_transform(t(Zq_corr), ref_umap_model)
    colnames(umap_query) = c("UMAP1", "UMAP2")
  }
  if (verbose)
    message("All done!")
  return(list(exp = exp_query, meta_data = metadata_query,
              Z = Zq_corr, Zq_pca = Z_pca_query, R = R_query, Xq = Xq,
              umap = umap_query))
}

# call clusters
callClusters  <- function(obj,
                          res=0.4,
                          k.near=30,
                          clustOB="svd",
                          cname="LouvainClusters",
                          min.reads=5e4,
                          m.clst=25,
                          threshold=5,
                          umap1="umap1",
                          umap2="umap2",
                          e.thresh=3,
                          cleanCluster=T,
                          cl.method=1,
                          svd_slotName="PCA",
                          umap_slotName="UMAP",
                          cluster_slotName="Clusters",
                          verbose=FALSE,
                          ...){
  
  # funcs
  filterSingle  <- function(pro,
                            k=50,
                            threshold=3,
                            type="umap",
                            m1="umap1",
                            m2="umap2"){
    
    # set column names
    if(type=="umap"){
      vars <- c(m1, m2)
    }else if(type=="pca"){
      vars <- colnames(pro)
    }
    
    # ensure that k is less than number of samples
    if(k > nrow(pro)){
      k <- nrow(pro)-1
    }
    
    # get nearest neighbors
    topk <- FNN::get.knn(pro[,vars], k=k)
    cell.dists <- as.matrix(topk$nn.dist)
    rownames(cell.dists) <- rownames(pro)
    colnames(cell.dists) <- paste0("k",seq(1:ncol(cell.dists)))
    aves <- apply(cell.dists, 1, mean)
    zscore <- as.numeric(scale(aves))
    names(zscore) <- rownames(pro)
    
    # thresholds
    p.zscore <- zscore[order(zscore, decreasing=T)]
    num.pass <- length(zscore[zscore < threshold])
    
    # filter
    prop.good <- zscore[zscore < threshold]
    ids <- names(prop.good)
    out <- pro[rownames(pro) %in% ids,]
    
    # return object
    return(out)
  }
  filtDistClst <- function(b,
                           umap1="umap1",
                           umap2="umap2",
                           threshold2=2){
    
    # iterate over each cluster
    clusts <- unique(b$seurat_clusters)
    out <- lapply(clusts, function(x){
      b.sub <- subset(b, b$seurat_clusters == x)
      b.umap <- b.sub[,c(umap1, umap2)]
      out.b.umap <- filterSingle(b.umap, k=25, threshold=threshold2, m1=umap1, m2=umap2)
      return(rownames(out.b.umap))
    })
    out <- do.call(c, out)
    b.out <- b[rownames(b) %in% as.character(out),]
    message("   * total number of cells surviving subcluster filtering = ", nrow(b.out))
    return(b.out)
    
  }
  
  # filter umap coordinates
  if(verbose){message(" - filtering outliers in UMAP manifold (z-score e.thresh = ", e.thresh, ") ...")}
  umap.original <- obj[[umap_slotName]]
  umap.filtered <- filterSingle(obj[[umap_slotName]], threshold=e.thresh, m1=umap1, m2=umap2)
  counts.filtered <- obj$counts[,rownames(umap.filtered)]
  meta.filtered <- obj$meta[rownames(umap.filtered),]
  pca.filtered <-  obj[[svd_slotName]][rownames(umap.filtered),]
  
  # run graph-based clustering
  if(verbose){message(" - creating seurat object for graph-based clustering ...")}
  
  # create Seurat object, find clusters
  sro <- SeuratObject::CreateSeuratObject(counts.filtered, min.cells=0, min.features=0)
  sro[["svd"]] <- CreateDimReducObject(embeddings = pca.filtered, key = "PC_", assay = DefaultAssay(sro))
  sro[["umap"]] <- CreateDimReducObject(embeddings=as.matrix(umap.filtered), key="UMAP_", assay=DefaultAssay(sro))
  sro <- AddMetaData(sro, meta.filtered)
  nn.eps.val <- 0
  n.starts <- 100
  
  sro <- FindNeighbors(sro, dims = 1:ncol(sro[[clustOB]]), reduction=clustOB, 
                       nn.eps=nn.eps.val, k.param=k.near, annoy.metric="euclidean")
  if(cl.method==4){
    graph.obj <- graph_from_adjacency_matrix(as(sro[[paste0(DefaultAssay(sro),"_snn")]], "dgCMatrix"),
                                             mode="undirected",
                                             weighted=T)
    
    graph.obj <- cluster_leiden(graph.obj, 
                                objective_function="modularity",
                                resolution_parameter=res)
    
    clustered <- membership(graph.obj)
    sro.meta <- data.frame(sro@meta.data)
    sro.meta$seurat_clusters <- factor(as.numeric(clustered[rownames(sro@meta.data)]))
    
  }else{
    sro <- FindClusters(sro, resolution=res, n.start=n.starts, algorithm=cl.method, ...)
    sro.meta <- data.frame(sro@meta.data)
    sro.meta$seurat_clusters <- factor(sro.meta$seurat_clusters)
  }  
  if(verbose){message(" - finished graph-based clustering ...")}
  
  # remove temp
  rm(umap.filtered)
  rm(counts.filtered)
  rm(meta.filtered)
  rm(pca.filtered)
  suppressMessages(gc())
  
  # remove outliers?
  if(cleanCluster){
    
    # verbose
    if(verbose){message("   * removing low quality clusters ...")}
    
    # prep
    sro.umap <- obj[[umap_slotName]][rownames(sro.meta),]
    colnames(sro.umap) <- c(umap1, umap2)
    sro.meta <- cbind(sro.meta, sro.umap)
    
    # filter by cluster size and number of cells
    agg.reads <- aggregate(sro.meta$nSites~sro.meta$seurat_clusters, FUN=sum)
    colnames(agg.reads) <- c("clusters","readDepth")
    clust.cnts <- table(sro.meta$seurat_clusters)
    agg.reads$num_cells <- clust.cnts[as.character(agg.reads$clusters)]
    agg.pass <- subset(agg.reads, agg.reads$num_cells >= m.clst & agg.reads$readDepth >= min.reads)
    sro.filt <- sro.meta[as.character(sro.meta$seurat_clusters) %in% as.character(agg.pass$clusters),]
    sro.filt$seurat_clusters <- droplevels(sro.filt$seurat_clusters)
    sro.filt$seurat_clusters <- as.numeric(factor(sro.filt$seurat_clusters))
    
    # remove outliers in the embedding
    if(verbose){message("   * filtering per-cluster outliers (z-score filtDistClst2 = ", threshold, ") ...")}
    sro.meta <- filtDistClst(sro.filt, umap1=umap1, umap2=umap2, threshold=threshold)
    sro.meta$seurat_clusters <- factor(sro.meta$seurat_clusters)
    
  }
  
  # filter by cluster size
  if(verbose){message(" - filtering clusters with low cell/read counts ...")}
  agg.reads <- aggregate(sro.meta$nSites~sro.meta$seurat_clusters, FUN=sum)
  colnames(agg.reads) <- c("clusters","readDepth")
  clust.cnts <- table(sro.meta$seurat_clusters)
  agg.reads$num_cells <- clust.cnts[as.character(agg.reads$clusters)]
  agg.pass <- subset(agg.reads, agg.reads$num_cells>=m.clst & agg.reads$readDepth>=min.reads)
  sro.filt <- sro.meta[sro.meta$seurat_clusters %in% agg.pass$clusters,]
  sro.filt$seurat_clusters <- droplevels(sro.filt$seurat_clusters)
  sro.filt$seurat_clusters <- as.numeric(as.factor(sro.filt$seurat_clusters))
  
  # rename output
  final.UMAP <- umap.original[rownames(sro.filt),]
  clusts <- sro.filt$seurat_clusters
  obj[[cluster_slotName]] <- obj$meta[rownames(sro.filt),]
  obj[[cluster_slotName]][,cname] <- factor(clusts)
  obj[[cluster_slotName]][,c(umap1)] <- final.UMAP[,c(umap1)]
  obj[[cluster_slotName]][,c(umap2)] <- final.UMAP[,c(umap2)]
  
  # return
  return(obj)
}

# save umap model
saveUMAP <- function(save_uwot_path, 
                     model, 
                     obj){
  
  if(file.exists(save_uwot_path)){
    message(paste('File already exists at that path... overwriting...'))
    file.remove(save_uwot_path)
  }
  
  # update object
  obj$umap$embedding <- model$embedding
  colnames(obj$umap$embedding) = c('UMAP1', 'UMAP2')
  model = uwot::save_uwot(model, file = save_uwot_path, unload = FALSE, verbose = FALSE)
  obj$save_uwot_path = save_uwot_path
  
  # verbose
  message(paste('Saved uwot model'))
  
  # return
  return(obj)
}


# load rds files and pre-processing --------------------------------------
out <- "maize_282.v8.3"
soc.obj <- readRDS("./QC_data/maize_282.merged.socObj.liberal.v8.rds")
meta.data <- read.table("./QC_data/All_pools.clusters.v1.tsv")
meta.data <- subset(meta.data, meta.data$type=="singlet" & (meta.data$gREF > meta.data$bREF) & meta.data$call==1 & meta.data$LLK_singlet > 3)
soc.obj$counts <- soc.obj$counts[,colnames(soc.obj$counts) %in% rownames(meta.data)]
soc.obj$meta <- soc.obj$meta[colnames(soc.obj$counts),]
meta.data <- meta.data[rownames(soc.obj$meta),]
soc.obj$meta$genotype <- meta.data$genotype
soc.obj$meta$library <- meta.data$library
soc.obj.raw <- soc.obj # save raw object


# get per cell feature counts --------------------------------------------
cell.counts <- log10(Matrix::colSums(soc.obj$counts))  # count number of features with Tn5 insertions per cell
cell.counts.z <- as.numeric(scale(cell.counts)) # convert features counts into Z-scores
cell.counts.threshold <- max(c((10^cell.counts[cell.counts.z < -1]), 1000)) # minimum feature counts (greater of 1 std or 1000)


# clean sparse counts matrix ---------------------------------------------
soc.obj <- cleanData(soc.obj, 
                     min.c=cell.counts.threshold,  # minimum number of accessible features per cell
                     min.t=0.0025, # minimum feature frequency across cells
                     max.t=0, # maximum feature frequency across cells
                     verbose=T)


# normalize with TFIDF ---------------------------------------------------
soc.obj <- tfidf(soc.obj, doL2=T)
number.sites <- ceiling(nrow(soc.obj$counts)*0.25)

# project with SVD -------------------------------------------------------
soc.obj <- reduceDims(soc.obj,
                      method="SVD", 
                      n.pcs=21,
                      cor.max=0.7,
                      num.var=number.sites,
                      verbose=T,
                      scaleVar=T,
                      doSTD=F,
                      doL1=F,
                      doL2=T,
                      refit_residuals=F)


# extract feature loadings and var/mean of tfidf
feat.data <- extractFeatureMetrics(soc.obj, num.var=number.sites)
ids <- rownames(soc.obj$PCA)


# remove batch effects with harmony --------------------------------------
ref.obj <- HarmonyMatrix(soc.obj$PCA, meta_data=soc.obj$meta, 
                         vars_use=c("library", "genotype"), 
                         do_pca=F,
                         theta=c(2, 2),
                         tau=c(5),
                         sigma=0.1,
                         lambda=c(0.1, 0.1),
                         nclust=50,
                         max.iter.cluster=100,
                         max.iter.harmony=30,
                         return_object=T)


# create compressed harmony reference
sym.ref <- symphony::buildReferenceFromHarmonyObj(ref.obj,
                                                  soc.obj$meta,
                                                  feat.data$features,
                                                  feat.data$loadings,               
                                                  verbose = TRUE,         
                                                  do_umap = F,       
                                                  umap_min_dist = 0.01,
                                                  save_uwot_path = paste0('./',out,'.uwot_model'))
soc.obj$PCA <- t(ref.obj$Z_corr)
soc.obj$PCA <- t(apply(soc.obj$PCA, 1, function(x){(x-mean(x))/sd(x)}))
colnames(soc.obj$PCA) <- paste0("PC_", 2:(ncol(soc.obj$PCA)+1))
rownames(soc.obj$PCA) <- ids


# reduce to 2-dimensions with UMAP ---------------------------------------
umap.modelout <- umap(soc.obj$PCA, 
                      metric="cosine",
                      n_neighbors=30,
                      a=2, b=0.75,
                      ret_model=T)


# add and save UMAP model
save_uwot_path = paste0('./',out,'.uwot_model')
sym.ref <- saveUMAP(save_uwot_path,
                    umap.modelout,
                    sym.ref)
saveRDS(sym.ref, file=paste0(out,".symphony_reference.rds"))


# update UMAP slot
soc.obj$UMAP <- umap.modelout$embedding
colnames(soc.obj$UMAP) <- c("umap1","umap2")
rownames(soc.obj$UMAP) <- rownames(soc.obj$PCA)


# identify clusters using neighborhood graph -----------------------------
soc.obj <- callClusters(soc.obj, 
                        res=2.0,
                        k.near=30,
                        verbose=T,
                        cleanCluster=F,
                        cl.method=4,
                        e.thresh=3,
                        threshold=3,
                        m.clst=50)


# plot cluster membership on UMAP embedding ------------------------------
pdf(paste0(out,".REF_CELLs.UMAP.clusters.pdf"), width=16, height=16)
plotUMAP(soc.obj, cluster_slotName="Clusters", cex=0.2)
dev.off()

pdf(paste0(out,".REF_CELLs.UMAP.library.pdf"), width=16, height=16)
plotUMAP(soc.obj, cluster_slotName="Clusters", column="library", cex=0.2)
dev.off()

pdf(paste0(out,".REF_CELLs.UMAP.nSites.pdf"), width=16, height=16)
plotUMAP(soc.obj, cluster_slotName="Clusters", column="log10nSites", cex=0.2)
dev.off()


# save data --------------------------------------------------------------
saveRDS(soc.obj, file=paste0(out,".REF_CELLs.processed.rds"))

# output text files
meta <- soc.obj$Clusters
rd <- soc.obj$PCA

# write data
write.table(meta, file=paste0(out, ".REF_CELLs.metadata.txt"), quote=F, row.names=T, col.names=T, sep="\t")
write.table(rd, file=paste0(out, ".REF_CELLs.reduced_dimensions.txt"), quote=F, row.names=T, col.names=T, sep="\t")


# map new data -----------------------------------------------------------
soc.obj.raw <- cleanData(soc.obj.raw,
                         min.c=100,  # minimum number of accessible features per cell
                         min.t=0, # minimum feature frequency across cells
                         max.t=0, # maximum feature frequency across cells
                         verbose=T)
new.obj <- cleanObj(soc.obj.raw, soc.obj)
new.obj <- projectTFIDF(new.obj, soc.obj, doL2=T)
splitter <- factor(sample(floor(nrow(new.obj$meta)/5000), nrow(new.obj$meta), replace=T))
ids <- split(rownames(new.obj$meta), splitter)
its <- 0
split.query <- lapply(ids, function(i){
  its <<- its + 1
  message(" - aligning query: ",names(ids)[its])
  message("   . query contrains ",nrow(new.obj$meta[i,]), " cells")
  mapQuery2(new.obj$residuals[,i],
            new.obj$meta[i,],
            sym.ref,
            svd_diag = soc.obj$PCA_model$d[soc.obj$PCA_model$keep_pcs],
            vars = c("library","genotype"),
            sigma = 0.01,
            do_umap = TRUE,
            scaleVar = F,
            doL2 = T,
            doSTD = F)

})

# merge mapped cells meta
num <- 0
all.meta <- lapply(split.query, function(x){
  num <<- num + 1
  colnames(x$umap) <- c("umap1_og","umap2_og")
  out <- cbind(x$meta,x$umap)
  out$map_type <- paste0("query",num)
  return(out)
})
all.meta <- do.call(rbind, all.meta)
rownames(all.meta) <- unlist(ids)
colnames(sym.ref$umap$embedding) <- c("umap1_og","umap2_og")
ref.meta <- cbind(sym.ref$meta, sym.ref$umap$embedding)
ref.meta$map_type <- "reference"
all.meta <- rbind(ref.meta, all.meta)

# merge mapped cells embeddings
all.embed <- lapply(split.query, function(x){
  t(x$Z)
})
all.embed <- do.call(rbind, all.embed)
all.embed <- rbind(t(sym.ref$Z_corr), all.embed)
colnames(all.embed) <- paste0("PC_", seq(1:ncol(all.embed)))

# make umap
all.umap <- all.meta[,c("umap1_og","umap2_og")]
rownames(all.umap) <- rownames(all.meta)
#all.meta[,c("umap1","umap2")] <- NULL

# merge counts
shared.sites <- intersect(rownames(new.obj$residuals), rownames(soc.obj$residuals))
res1 <- new.obj$residuals[shared.sites,]
res2 <- soc.obj$residuals[shared.sites,]
counts <- cbind(res1, res2)
shared <- intersect(rownames(all.meta), colnames(counts))
counts <- counts[,shared]
all.embed <- all.embed[shared,]
all.meta <- all.meta[shared,]
all.umap <- all.umap[shared,]

# new object
m.obj <- list(PCA=as.matrix(all.embed),
              meta=all.meta,
              UMAP=all.umap,
              counts=counts)

# L2 normalize reduced embeddings
m.obj$PCA.l2 <- m.obj$PCA #t(apply(m.obj$PCA, 1, function(x){(x-mean(x,na.rm=T))/sd(x, na.rm=T)}))
m.obj$UMAP.l2 <- umap(m.obj$PCA.l2,
                      metric="correlation",
                      n_neighbors=15,
                      a=1.9, b=0.6)
colnames(m.obj$UMAP.l2) <- c("umap1","umap2")
rownames(m.obj$UMAP.l2) <- rownames(m.obj$PCA.l2)


# re-identify clusters ----------------------------------------------------
m.obj$PCA.l2 <- m.obj$PCA #t(apply(m.obj$PCA, 1, function(x){(x-mean(x))/sd(x)}))
m.obj <- callClusters(m.obj,
                      res=1.5,
                      k.near=10,
                      verbose=T,
                      cleanCluster=F,
                      cl.method=4,
                      e.thresh=5,
                      threshold=5,
                      m.clst=50,
                      svd_slotName="PCA.l2",
                      umap_slotName="UMAP.l2")


# plot cluster membership on UMAP embedding ------------------------------
dims <- 16
pdf(paste0(out,".ALL_CELLs.UMAP.clusters.pdf"), width=dims, height=dims)
plotUMAP(m.obj, cluster_slotName="Clusters", cex=0.2)
dev.off()

pdf(paste0(out,".ALL_CELLs.UMAP.library.pdf"), width=dims, height=dims)
plotUMAP(m.obj, cluster_slotName="Clusters", column="library", cex=0.2)
dev.off()

pdf(paste0(out,".ALL_CELLs.UMAP.nSites.pdf"), width=dims, height=dims)
plotUMAP(m.obj, cluster_slotName="Clusters", column="log10nSites", cex=0.2)
dev.off()


# write data -------------------------------------------------------------
write.table(m.obj$Clusters, file=paste0(out, ".ALL_CELLs.metadata.txt"), quote=F, row.names=T, col.names=T, sep="\t")
write.table(m.obj$PCA.l2[rownames(m.obj$Clusters),], file=paste0(out, ".ALL_CELLs.reduced_dimensions.txt"), quote=F, row.names=T, col.names=T, sep="\t")
saveRDS(split.query, file=paste0(out,".ALL_CELLs.aligned_query.rds"))
saveRDS(m.obj, file=paste0(out,".ALL_CELLs.processed.rds"))
