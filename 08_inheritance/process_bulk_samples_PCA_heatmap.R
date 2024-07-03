## plot PCA/heatmap of bulk ATAC samples ##

# load libraries
library(viridis)
library(gplots)
library(RColorBrewer)
library(wesanderson)

# load data
a <- read.table("filtered_samples_CPM_in_peaks.tsv", header=T)

# mats
z <- as.matrix(t(scale(t(a))))
z <- z[complete.cases(z),]
pca <- prcomp(t(z), center=F, scale.=F)
df <- as.data.frame(do.call(rbind, strsplit(colnames(z),"_")))
rownames(df) <- paste(df$V1,df$V2,sep="_")
colnames(z) <- rownames(df)
df$oo <- factor(df$V1, levels=c("B73", "Ki3", "Oh43", "BxK", "BxO", "OxK"))
df <- df[order(df$oo, decreasing=F),]
z <- z[,rownames(df)]

# plot info
cols1 <- wes_palette("Darjeeling1", 5)
cols2 <- wes_palette("Darjeeling2", 5)
cols <- c(cols1, cols2[1:2])
type <- unique(df$V1)
type <- type[c(1,4,5,2,3,6)]
names(cols) <- type
df$col <- cols[df$V1]

# plot
pdf("F1_P.PCA.pdf", width=5, height=5)
plot(pca$x[,1:2], col=df$col, pch=16, xlab="PC1", ylab="PC2")
grid(lty=1)
legend("topright", legend=type, fill=cols, border=NA)
dev.off()

# process z-score matrix
mods <- lapply(seq(1:nrow(z)), function(x){
  if((x %% 1000)==0){message(" - iterated over ",x," records...")}
  ddf <- cbind(df, as.numeric(z[x,]))
  colnames(ddf)[ncol(ddf)] <- c("accessibility")
  colnames(ddf)[1] <- c('genotype')
  mod <- lm(accessibility~genotype, data=ddf)
  res <- anova(mod)
  pval <- res$`Pr(>F)`[1]
  return(pval)
})
mods <- unlist(mods)
score <- data.frame(pval=mods, row.names=rownames(z))
score$fdr <- p.adjust(score$pval, method="fdr")
score$keeper <- ifelse(score$fdr < 0.05, 1, 0)

# heatmap
cols <- wes_palette(name="Zissou1", 100, type="continuous")
col.col <- as.character(cols)

# ref vals
zz <- z[score$keeper==1,]
zz[zz > 3] <- 3
zz[zz < -3] <- -3

# plot
pdf("ACR_x_sample_heatmap.FDR05.pdf", width=8, height=8)
heatmap.2(zz, 
          Colv=F, Rowv=T, 
          labRow=F,
          trace='none', 
          useRaster=T, col=col.col, dendrogram='row')
dev.off()