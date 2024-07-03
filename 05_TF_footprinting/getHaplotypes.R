library(impute)
library(mclust)

# load data
a <- readRDS("caQTL-ACRs.celltypes.SNVs.imputed.rds")

# iterate over each ACR
doCluster <- T
acrs <- unique(a$acrID)
its <- 0
outs <- lapply(acrs, function(x){
  
  # verbose
  its <<- its + 1
  if((its %% 100)==0){message(" - iterated over ", its, " ACRs...")}
  
  # subset by ACR
  df <- subset(a, a$acrID==x)
  z <- as.matrix(df[,5:ncol(df)])
  
  # 2 or more snps
  if(nrow(z)>1){
    
    # mask completely missing individuals
    all.ids <- colnames(z)
    keep.ids <- apply(z, 2, function(y){
      if(sum(is.na(y))>(nrow(z)*0.95)){
        return(F)
      }else{
        T
      }
    })
    z <- z[,keep.ids]
    
    # estimate distance and cluster
    cors <- as.matrix(dist(t(z)))
    km <- suppressWarnings(suppressMessages(Mclust(cors, G=1:10, verbose=F)))
    clusters <- km$classification
    counts <- table(clusters)
    counts <- counts[order(counts, decreasing=T)]
    if(length(counts) > 2){
      nottop.2 <- names(counts)[3:length(counts)]
      clusters[clusters %in% nottop.2] <- NA
    }
    if(!identical(names(counts)[1:2], c("1","2"))){
      two.ids <- names(clusters)[which(clusters == names(counts)[2])]
      one.ids <- names(clusters)[which(clusters == names(counts)[1])]
      clusters <- ifelse(names(clusters) %in% one.ids, 1,
                         ifelse(names(clusters) %in% two.ids, 2, clusters))
      names(clusters) <- names(km$classification)
    }
    clusters <- clusters[all.ids]
    names(clusters) <- all.ids
    
    # SNV profile cluster 1
    hap1 <- rowMeans(z[,!is.na(names(clusters[clusters==1]))])
    hap2 <- rowMeans(z[,!is.na(names(clusters[clusters==2]))])
    ddf <- df[,1:4]
    ddf$hap1 <- hap1
    ddf$hap2 <- hap2
    
    # cluster SNVs
    non.na <- clusters[!is.na(clusters)]
    zz <- z[,names(non.na)[order(non.na)]]
    
    # return data
    if(doCluster){
      return(clusters)
    }else{
      return(ddf)
    }
    
  }else{
    
    # cluster SNP
    all.ids <- colnames(z)
    km <- suppressWarnings(suppressMessages(Mclust(z, G=1:10, verbose=F)))
    clusters <- km$classification
    counts <- table(clusters)
    counts <- counts[order(counts, decreasing=T)]
    if(length(counts) > 2){
      nottop.2 <- names(counts)[3:length(counts)]
      clusters[clusters %in% nottop.2] <- NA
    }
    if(!identical(names(counts)[1:2], c("1","2"))){
      two.ids <- names(clusters)[which(clusters == names(counts)[2])]
      one.ids <- names(clusters)[which(clusters == names(counts)[1])]
      clusters <- ifelse(names(clusters) %in% one.ids, 1,
                         ifelse(names(clusters) %in% two.ids, 2, clusters))
      names(clusters) <- names(km$classification)
    }
    clusters <- clusters[all.ids]
    names(clusters) <- all.ids
    
    # SNV profile cluster 1
    hap1 <- mean(z[!is.na(names(clusters[clusters==1]))])
    hap2 <- mean(z[!is.na(names(clusters[clusters==2]))])
    ddf <- df[,1:4]
    ddf$hap1 <- hap1
    ddf$hap2 <- hap2
    
    # return data
    if(doCluster){
      return(clusters)
    }else{
      return(ddf)
    }
    
  }
  
})
outs <- do.call(rbind, outs)