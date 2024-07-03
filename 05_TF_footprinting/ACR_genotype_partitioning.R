## cluster ACR SNV genotypes ##

# load libraries
library(impute)
library(mclust)

# load data
a <- read.table(gzfile("caQTL-ACRs.celltypes.SNVs.bed.gz"), header=T, check.names=F)

# update ACR id
a$QUAL <- paste(a$acrChr, (a$acrStr+1000), (a$acrEnd-1000), sep="_")
a$ID <- paste(a$CHROM, a$POS, sep="_")
a[,1:5] <- NULL
a$FILTER <- NULL
a$INFO <- NULL
a$FORMAT <- NULL
colnames(a)[1] <- "snpID"
colnames(a)[4] <- "acrID"
a <- a[,c(1,4,2,3,5:ncol(a))]

# simplify genotype
b <- apply(a[,5:ncol(a)], 2, function(x){
  xx <- ifelse(grepl('0/0', x), 0,
               ifelse(grepl('0/1', x), 1,
                      ifelse(grepl('1/1', x), 2, NA)))
  return(xx)
})

# remove highly heterozygous sites
het.freq <- t(apply(b, 1, function(x){
  fac <- factor(x, levels=c(0,1,2))
  props <- prop.table(table(fac))
  return(props)
}))
keepers <- het.freq[,2] < 0.1
a <- a[keepers,1:4]
b <- b[keepers,]

# impute missing data
imp <- suppressWarnings(suppressMessages(impute.knn(b, k=10, rowmax=0.99, colmax=0.99)$data))
a <- cbind(a, imp)

# iterate over ACRs and select top two haplotypes
acrs <- unique(a$acrID)
its <- 0
doCluster <- F
outs <- lapply(acrs, function(x){
  
  # which output
  doCluster <- F
  
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

# save
if(doCluster){
  
  # add rownames
  rownames(outs) <- acrs
  
  # clean up tags
  new <- t(apply(outs, 1, function(x){
    counts <- table(x)
    ids <- names(counts)
    if(!identical(ids, c("1", "2"))){
      id1 <- ids[1]
      xx <- ifelse(x==id1, "1", "2")
      return(xx)
    }else{
      return(x)
    }
  }))
  new <- data.matrix(new)
  write.table(new, file="ACR_haplotype_groups.txt", quote=F, row.names=T, col.names=T, sep="\t")
  
}else{
  
  # output SNV data
  write.table(outs, file="SNV_haplotype_groups.txt", quote=F, row.names=T, col.names=T, sep="\t")
  
}