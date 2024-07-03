# set eFDR for coACRs

# libraries
library(scales)
library(png)

# load data
args <- commandArgs(T)
if(length(args) != 3){stop("Rscript eFDR_coACRs.R <coaccess.rds> <perm_coaccess.rds> <prefix>")}
obsD <- as.character(args[1])
permD <- as.character(args[2])
prefix <- as.character(args[3])

# function
plotData <- function(obs){
    
    # plot
    obs <- obs[order(obs$call, decreasing=T),]
    png("All_coACR_correlation_v_length.png", width=6, height=5.5, unit="in", res=300)
    plot(log10(obs$distance[obs$coaccess.adj > 0]),
         obs$coaccess.adj[obs$coaccess.adj > 0], 
         col=ifelse(obs$call[obs$coaccess.adj > 0]==1, "dodgerblue4", "grey80"), 
         cex=0.2, pch=16,
         xlab="coACR length (log10)",
         ylab="Adjusted co-accessibility")
    dev.off()
    
    pass <- subset(obs, obs$call==1)
    fail <- subset(obs, obs$call==0)
    
    den1 <- density(log10(pass$distance))
    den2 <- density(log10(fail$distance))
    
    pdf("All_coACR_length_dist.pdf", width=6, height=5.5)
    plot(den1, col="dodgerblue4", ylim=c(0, max(den1$y, den2$y)), lwd=2)
    lines(den2, col="grey80", lwd=2)
    dev.off()
    
}

# load data
a <- readRDS(obsD)
b <- readRDS(permD)
#acrs <- read.table("/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/ACR_analysis/consistentACRs.distTSS.bed")

# filter on eFDR 
efdr <- 0.05
if(class(a)=="data.frame"){
    
    # load data
    cluster <- b
    obs <- a
    
    # rename columns
    colnames(cluster)[3] <- "coaccess"
    colnames(obs)[3] <- "coaccess"
    
    # select the same coACRs
    rownames(obs) <- paste(obs$Var1, obs$Var2,sep="-")
    rownames(cluster) <- paste(cluster$Var1, cluster$Var2,sep="-")
    shared <- intersect(rownames(obs), rownames(cluster))
    obs <- obs[shared,]
    cluster <- cluster[shared,]
    rownames(obs) <- seq(1:nrow(obs))
    rownames(cluster) <- seq(1:nrow(cluster))
    
    # estimate coACR distance
    df1 <- data.frame(do.call(rbind, strsplit(as.character(obs$Var1),"_")))
    df2 <- data.frame(do.call(rbind, strsplit(as.character(obs$Var2),"_")))
    df1$X2 <- as.numeric(as.character(df1$X2))
    df1$X3 <- as.numeric(as.character(df1$X3))
    df2$X2 <- as.numeric(as.character(df2$X2))
    df2$X3 <- as.numeric(as.character(df2$X3))
    df1$ave <- as.integer((df1$X2+df1$X3)/2)
    df2$ave <- as.integer((df2$X2+df2$X3)/2)
    obs$distance <- abs(df2$ave-df1$ave)
    
    message(" - Num OBS = ", nrow(obs), " | Num EXP = ", nrow(cluster))
    
    # penalize coaccessibility by distance
    obs$coaccess.adj <- obs$coaccess * (sqrt(2000)/sqrt(obs$distance))
    cluster$coaccess.adj <- cluster$coaccess * (sqrt(2000)/sqrt(obs$distance))
    
    # subset cluster
    pos <- subset(cluster, cluster$coaccess.adj >= 0)
    pos.threshold <- signif(quantile(pos$coaccess.adj, 1-efdr))
    message(" - threshold (+) = ", pos.threshold)#," | (-) = ", neg.threshold)
    obs$adj.residual <- (obs$coaccess.adj - cluster$coaccess.adj)
    obs$coaccess.residual <- obs$coaccess - cluster$coaccess
    obs$call <- ifelse(obs$coaccess.adj >= pos.threshold & obs$adj.residual > 0, 1, 0)# | obs$value <= neg.threshold)
    ratio <- mean(obs$call)
    ratio <- signif(ratio, digits=6)
    message(" - filtered coACR set = ", sum(obs$call), " | ", ratio)
    message("")
    
    # process coACRs
    obs$coACR <- paste(obs$Var1, obs$Var2,sep="-")
    obs$Var1 <- NULL
    obs$Var2 <- NULL
    obs <- obs[,c(6,1,3,4,2,5)]
    filtered <- obs #subset(obs, obs$coaccess > 0)
    
}else{
    ids <- names(a)

    filtered <- lapply(ids, function(x){
        
        # load data
        cluster <- b[[x]]
        obs <- a[[x]]
        
        # rename columns
        colnames(cluster)[3] <- "coaccess"
        colnames(obs)[3] <- "coaccess"
        
        # select the same coACRs
        rownames(obs) <- paste(obs$Var1, obs$Var2,sep="-")
        rownames(cluster) <- paste(cluster$Var1, cluster$Var2,sep="-")
        shared <- intersect(rownames(obs), rownames(cluster))
        obs <- obs[shared,]
        cluster <- cluster[shared,]
        rownames(obs) <- seq(1:nrow(obs))
        rownames(cluster) <- seq(1:nrow(cluster))
        
        # estimate coACR distance
        df1 <- data.frame(do.call(rbind, strsplit(as.character(obs$Var1),"_")))
        df2 <- data.frame(do.call(rbind, strsplit(as.character(obs$Var2),"_")))
        df1$X2 <- as.numeric(as.character(df1$X2))
        df1$X3 <- as.numeric(as.character(df1$X3))
        df2$X2 <- as.numeric(as.character(df2$X2))
        df2$X3 <- as.numeric(as.character(df2$X3))
        df1$ave <- as.integer((df1$X2+df1$X3)/2)
        df2$ave <- as.integer((df2$X2+df2$X3)/2)
        obs$distance <- abs(df2$ave-df1$ave)
        
        message(" - Num OBS = ", nrow(obs), " | Num EXP = ", nrow(cluster))
        
        # penalize coaccessibility by distance
        obs$coaccess.adj <- obs$coaccess * (sqrt(2000)/sqrt(obs$distance))
        cluster$coaccess.adj <- cluster$coaccess * (sqrt(2000)/sqrt(obs$distance))
        
        # subset cluster
        pos <- subset(cluster, cluster$coaccess.adj >= 0)
        pos.threshold <- signif(quantile(pos$coaccess.adj, 1-efdr))
        message(" - thresholds for cluster ",x," | (+) = ", pos.threshold)#," | (-) = ", neg.threshold)
        obs$adj.residual <- (obs$coaccess.adj - cluster$coaccess.adj)
	obs$coaccess.residual <- obs$coaccess - cluster$coaccess
        obs$call <- ifelse(obs$coaccess.adj >= pos.threshold & obs$adj.residual > 0, 1, 0)# | obs$value <= neg.threshold)
        ratio <- mean(obs$call, na.rm=T)
        ratio <- signif(ratio, digits=6)
        message(" - filtered coACR set = ", sum(obs$call, na.rm=T), " | ", ratio)
        message("")
        
        # process coACRs
        obs$coACR <- paste(obs$Var1, obs$Var2,sep="-")
        obs$Var1 <- NULL
        obs$Var2 <- NULL
        obs <- obs[,c(8,1,4,5,3,2,6,7)]
        #obs <- subset(obs, obs$coaccess > 0)
        return(obs)
    })
    filtered <- do.call(rbind, filtered)
}

# plot
#plotData(filtered)

# anchor to tss
#tss <- subset(acrs, acrs$V4 < 50)
#tss <- tss[!duplicated(tss[,1:3]),]
#rownames(tss) <- paste(tss$V1,tss$V2,tss$V3,sep="_")

# check cluster 1
#obs <- subset(filtered, filtered$LouvainClusters==1 & filtered$call==1)
#obs <- cbind(obs, data.frame(do.call(rbind, strsplit(as.character(obs$coACR),"-"))))
#colnames(obs)[(ncol(obs)-1):ncol(obs)] <- c("acrID1", "acrID2")
#obs$acrID1 <- as.character(obs$acrID1)
#obs$acrID2 <- as.character(obs$acrID2)
#tss.obs <- obs[obs$acrID1 %in% rownames(tss) | obs$acrID2 %in% rownames(tss),]

# plot
#den1 <- density(log10(acrs$V4+1))
#den2 <- density(log10(tss.obs$distance))

#pdf("TSS_anchored.vline50kb.pdf", width=6, height=6)
#plot(den1, ylim=c(0, max(den1$y, den2$y)), col="grey75", xlim=c(1, 6))
#lines(den2, col="dodgerblue4")
#abline(v=log10(50000))
#grid()
#dev.off()

# final set
#final <- subset(filtered, filtered$call==1)
write.table(filtered, file=paste0(prefix,".coACRs.raw.txt"), quote=F, row.names=T, col.names=T, sep="\t")
