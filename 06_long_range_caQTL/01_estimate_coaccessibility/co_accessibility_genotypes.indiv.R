## co-accessibility ##

# load libraries
library(Matrix)
library(FNN)
library(gtools)
library(reshape2)
library(edgeR)
library(parallel)

# functions
coAccessibility <- function(obj, 
                            k=NULL, 
                            n.seed.cells=NULL, 
                            clusterID="celltype_full", 
                            embedding="PCA", 
                            threads=1, 
                            max.dist=5e5, 
                            split_data=F){
    
    # order peaks
    obj$peak_counts <- obj$peak_counts[mixedorder(rownames(obj$peak_counts), decreasing=F),]
    peaks <- data.frame(do.call(rbind, strsplit(rownames(obj$peak_counts),"_")))
    peaks$id <- rownames(obj$peak_counts)
    peaks$X1 <- as.character(peaks$X1)
    peaks$X2 <- as.numeric(as.character(peaks$X2))
    peaks$X3 <- as.numeric(as.character(peaks$X3))
    colnames(peaks) <- c("chr","str","end","id")
    peaks$ave <- as.integer((peaks$str+peaks$end)/2)
    peaks <- peaks[order(peaks$chr, peaks$str, peaks$end, decreasing=F),]
    rownames(peaks) <- peaks$id
    
    # iterate over clusters
    if(split_data){
        outs <- mclapply(unique(obj[["Clusters"]][,clusterID]), function(x){
            
            # verbose
            message(" - creating pseudocells for cluster ", x)
            
            # select cells in cluster
            df <- subset(obj[["Clusters"]], obj[["Clusters"]][,clusterID] == x)
            ids <- rownames(df)
            cnts <- obj$peak_counts[,ids]
            cnts <- cnts[Matrix::rowSums(cnts)>0,]
            shared <- intersect(rownames(cnts), rownames(peaks))
            peaks.s <- peaks[shared,]
            p.ids <- rownames(peaks.s)
            
            # iterated
            message(" - calculating co-accessibility across nuclei for cluster ", x)
            co.acrs <- lapply(seq(1:(nrow(peaks.s)-1)), function(y){
                if((y %% 1000) == 0){message(' - iterated over ', y,' ACRs | cluster ', x)}
                sites <- c(p.ids[y])
                for (i in (y+1):nrow(peaks.s)){
                    if(abs((peaks.s$ave[i] - peaks.s$ave[y]) < max.dist) & (peaks.s$chr[i] == peaks.s$chr[y])){ 
                        if((abs(peaks.s$ave[i] - peaks.s$ave[y]) > 2000)){
                            sites <- c(sites, p.ids[i])
                        }else{
                            next
                        }
                    }else{
                        break
                    }
                }
                keepers <- paste(rep(p.ids[y], length(sites)-1),
                                 sites[2:length(sites)],
                                 sep="-")
                if(length(sites) > 1){
                    val <- cor(as.matrix(t(cnts[sites,])), method="spearman")
                    val[lower.tri(val, diag = TRUE)] <- NA
                    val <- melt(val, na.rm=T)
                    rownames(val) <- paste(val$Var1, val$Var2, sep="-")
                    val <- val[keepers,]
                    rownames(val) <- seq(1:nrow(val))
                    return(val)
                }
            })
            message(" - aggregating co-accessibility scores for cluster ", x)
            co.acrs <- do.call(rbind, co.acrs)
            co.acrs <- co.acrs[order(co.acrs$Var1, co.acrs$Var2, decreasing=F),]
            co.acrs <- co.acrs[!duplicated(co.acrs),]
            co.acrs[,clusterID] <- x
            #saveRDS(co.acrs, file=paste0(x,"_coACRs.rds"))
            return(co.acrs)
        }, mc.cores=threads)
        
        # return output
        names(outs) <- unique(obj[["Clusters"]][,clusterID])
        return(outs)
        
    }else{
        
        # make pseudocells
        outs <- mclapply(unique(obj[["Clusters"]][,clusterID]), function(x){
            
            # verbose
            message(" - creating pseudocells for cluster ", x)
            
            # select cells in cluster
            df <- subset(obj[["Clusters"]], obj[["Clusters"]][,clusterID] == x)
            ids <- rownames(df)
            embed <- obj[[embedding]][ids,]
            if(is.null(k)){
                k <- as.integer(sqrt(nrow(df)))
            }
            if(is.null(n.seed.cells)){
                n.seed.cells <- as.integer((length(ids)/(0.9*k)))
            }
            message(" - k = ", k, " | n.seeds = ",n.seed.cells," | cluster ", x)
            
            # select territory seed cells
            s.cells <- sample(ids, n.seed.cells)
            
            # get knn
            knn.out <- FNN::get.knn(embed, k=k)$nn.index
            rownames(knn.out) <- rownames(embed)
            colnames(knn.out) <- paste0("KNN", seq(1:ncol(knn.out)))
            knn.out <- apply(knn.out, 2, function(z){
                rownames(knn.out)[z]
            })
            rownames(knn.out) <- rownames(embed)
            
            # make pseudocells
            s.knn <- knn.out[s.cells,]
            pcells <- lapply(seq(1:nrow(s.knn)), function(z){
                Matrix::rowSums(obj$peak_counts[,c(as.character(s.knn[z,]),rownames(s.knn)[z])])
            })
            pcells <- do.call(cbind, pcells)
            colnames(pcells) <- paste0("pseudo_",x,".",seq(1:ncol(pcells)))
            pcells <- cpm(pcells, log=F)
            return(pcells)
        })
        pcells <- do.call(cbind, outs)
        
        # run correlations
        shared.ids <- intersect(rownames(peaks), rownames(pcells))
        pcells <- pcells[shared.ids,]
        peaks.s <- peaks[shared.ids,]
        p.ids <- rownames(peaks.s)
        message(" - calculating co-accessibility across pseudocells")
        co.acrs <- lapply(seq(1:(nrow(peaks.s)-1)), function(y){
            if((y %% 1000) == 0){message(' - iterated over ', y,' ACRs')}
            sites <- c(p.ids[y])
            for (i in (y+1):nrow(peaks.s)){
                if(abs((peaks.s$ave[i] - peaks.s$ave[y]) < max.dist) & (peaks.s$chr[i] == peaks.s$chr[y])){ 
                    if((abs(peaks.s$ave[i] - peaks.s$ave[y]) > 2000)){
                        sites <- c(sites, p.ids[i])
                    }else{
                        next
                    }
                }else{
                    break
                }
            }
            keepers <- paste(rep(p.ids[y], length(sites)-1),
                             sites[2:length(sites)],
                             sep="-")
            if(length(sites) > 1){
                val <- cor(as.matrix(t(pcells[sites,])), method="spearman")
                val[lower.tri(val, diag = TRUE)] <- NA
                val <- melt(val, na.rm=T)
                rownames(val) <- paste(val$Var1, val$Var2, sep="-")
                val <- val[keepers,]
                rownames(val) <- seq(1:nrow(val))
                return(val)
            }
        })
        message(" - aggregating co-accessibility scores")
        co.acrs <- do.call(rbind, co.acrs)
        co.acrs <- co.acrs[order(co.acrs$Var1, co.acrs$Var2, decreasing=F),]
        co.acrs <- co.acrs[!duplicated(co.acrs),]
        return(co.acrs)
        
    }
}

# load data
message(" - loading input data ...")
dir <- "/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/step0_QC_data/"
acrdir <- "/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/sparse/acrs/"
metadata <- read.table(paste0(dir,"maize_282.v8.3.ALL_CELLs.metadata.leiden.pruned.subclustered.txt"))
peaks <- readRDS(paste0(acrdir,"mergedPools.celltype_genotype.ACRs.109k.sparse.rds"))
pcs <- read.table(paste0(dir,"maize_282.v8.3.ALL_CELLs.reduced_dimensions.txt"))

# process
soc.obj <- list(peak_counts=peaks,
                Clusters=metadata,
                PCA=pcs)

# clean
rm(peaks)
rm(pcs)
rm(metadata)

# run
coaccess <- coAccessibility(soc.obj, clusterID="genotype", threads=10, split_data=T)
saveRDS(coaccess, file="maize_282.co_accessibility.per_nucleus_genotypes.indiv.rds")

# permute
c.ids <- sample(colnames(soc.obj$peak_counts))
r.ids <- sample(rownames(soc.obj$peak_counts))
colnames(soc.obj$peak_counts) <- c.ids
rownames(soc.obj$peak_counts) <- r.ids

# run
coaccess <- coAccessibility(soc.obj, clusterID="genotype", threads=10, split_data=T)
saveRDS(coaccess, file="maize_282.co_accessibility.permuted.per_nucleus_genotypes.indiv.rds")
