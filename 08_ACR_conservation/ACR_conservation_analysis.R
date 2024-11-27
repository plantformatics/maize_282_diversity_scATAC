## id PAV acrs ##

# load libraries
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(gtools)
library(edgeR)
library(zFPKM)
library(MASS)
library(scales)
library(uwot)
library(reshape)
library(png)
library(Matrix)
library(glmnet)
library(parallel)

# functions
plot_cbd <- function(x1,x2,
                     ylim=c(min(x2),max(x2)),
                     xlim=c(min(x1),max(x1)),
                     xlab="",ylab="",main="",
                     cex=0.5, fit=F, bwf=5, nbin=300,
                     colP=NULL, rasterize=F){
    
    .normalize <- function(x){
        (x - min(x))/(max(x)-min(x))
    }
    if(is.null(colP)){
        
        colP <- colorRampPalette(c("grey75","grey70","darkorchid4","firebrick3","darkorange",
                           "gold1","yellow"), bias=1.5)(256)
        
    }else{
        colP <- colP(256)
    }
    bww <- (max(x1)-min(x1))/bwf
    bwh <- (max(x2)-min(x2))/bwf
    df <- data.frame(x1,x2)
    x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")),
                  nbin=nbin, bandwidth=c(bww, bwh))
    df$dens <- col2rgb(x)[1,] + 1L
    cols <- colP
    df$col <- cols[df$dens]
    df$opa <- .normalize(df$dens)
    
    # if raster
    if(rasterize){
        plot(x2~x1, data=df[order(df$dens),], 
             ylim=ylim,xlim=xlim,pch=16,col=alpha(col,1),
             cex=cex,xlab=xlab,ylab=ylab,
             main=main, type="n")
        coords <- par("usr")
        gx <- grconvertX(coords[1:2], "user", "inches")
        gy <- grconvertY(coords[3:4], "user", "inches")
        width <- max(gx) - min(gx)
        height <- max(gy) - min(gy)
        tmp <- tempfile()
        png(tmp, width = width, height = height, units = "in", res = 300, bg = "transparent")
        par(mar=c(0,0,0,0))
        plot.new()
        plot.window(coords[1:2], coords[3:4], mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
        points(x2~x1, data=df[order(df$dens),], 
               pch=16,col=alpha(col,1),
               cex=cex)
        dev.off()
        panel <- readPNG(tmp)
        rasterImage(panel, coords[1], coords[3], coords[2], coords[4])
        
    }else{
        plot(x2~x1, data=df[order(df$dens),], 
             ylim=ylim,xlim=xlim,pch=16,col=alpha(col,1),
             cex=cex,xlab=xlab,ylab=ylab,
             main=main)
    }
    grid(lty=1)
    
    # fit regression?
    if(fit){
        
        mod <- lm(x2~x1)
        abline(mod, col="firebrick3", lwd=2)
        
    }
    
    cors <- cor(x2,x1)
    mtext(paste0("PCC = ", signif(cors, digits=3)))
}
plot_raster <- function(x, y,
                        ylim=c(min(y),max(y)),
                        xlim=c(min(1),max(x)),
                        xlab="",ylab="",main="",
                        cex=0.5, pch=16,
                        col="black"){
    
    # plot first 
    plot(y~x, 
         ylim=ylim,xlim=xlim,
         pch=pch,col=col,
         cex=cex,xlab=xlab,ylab=ylab,
         main=main, type="n")
    coords <- par("usr")
    gx <- grconvertX(coords[1:2], "user", "inches")
    gy <- grconvertY(coords[3:4], "user", "inches")
    width <- max(gx) - min(gx)
    height <- max(gy) - min(gy)
    tmp <- tempfile()
    png(tmp, width = width, height = height, units = "in", res = 300, bg = "transparent")
    par(mar=c(0,0,0,0))
    plot.new()
    plot.window(coords[1:2], coords[3:4], mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
    points(y~x, 
           pch=pch,col=col,
           cex=cex)
    dev.off()
    panel <- readPNG(tmp)
    rasterImage(panel, coords[1], coords[3], coords[2], coords[4])
}
calcOverlapACR <- function(x){
    
    s1 <- x$V2
    e1 <- x$V3
    s2 <- x$V8
    e2 <- x$V9
    len <- ifelse(s1 < s2 & s2 < e1, e1-s2,
                  ifelse(s1 > s2 & e1 < e2, e1-s1,
                         ifelse(s2 < s1 & s1 < e2 & e1 > e2, e2-s1, e1-s1)))
    acr.len <- e1-s1
    frac <- len/acr.len
    return(frac)
    
}

# load data
scatac <- readRDS("/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/celltype_genotype_bulks/celltype_genotype.counts.raw.rds")
wgs <- read.table("all_genotype_celltype.WGS_counts.raw.txt")
cov <- read.table("all_genotype_celltype.WGS_counts.raw_cov.txt")
inbred.data <- read.table("Inbred_pedigree.v2.txt")
acr.data <- read.table("../celltype_genotype.ACRs.109k.features.STARR.bed")
phylo.data <- read.table("/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/step15_pop_gen_stats/Phylo_ACR_weighted_mean.txt")
teo <- read.table("Teosinte_ACR_coverage.bed")

# add tropical/nontropical column
inbred.data$type <- ifelse(inbred.data$TS > 0.9, "tropical", ifelse(inbred.data$TS < 0.1, "nontropical", "mixed"))
rownames(inbred.data) <- gsub("Goodman-Buckler","",rownames(inbred.data))

# reformat scATAC colnames
#colnames(scatac) <- paste0("282set_",colnames(scatac))
colnames(scatac) <- gsub("Goodman.Buckler","",colnames(scatac))
colnames(scatac) <- gsub("33\\.16","33-16",colnames(scatac))
colnames(scatac) <- gsub("38\\.11","38-11",colnames(scatac))
colnames(scatac) <- gsub("A441\\.5","A441-5",colnames(scatac))
colnames(scatac) <- gsub("CH701\\.30","CH701-30",colnames(scatac))
colnames(scatac) <- gsub("CI187\\.2","CI187-2",colnames(scatac))

# reformat wgs colnames
colnames(wgs) <- gsub("X282set_","282set_",colnames(wgs))
colnames(wgs) <- gsub("Goodman.Buckler","",colnames(wgs))
colnames(wgs) <- gsub("33\\.16","33-16",colnames(wgs))
colnames(wgs) <- gsub("38\\.11","38-11",colnames(wgs))
colnames(wgs) <- gsub("A441\\.5","A441-5",colnames(wgs))
colnames(wgs) <- gsub("CH701\\.30","CH701-30",colnames(wgs))
colnames(wgs) <- gsub("CI187\\.2","CI187-2",colnames(wgs))

# reformat cov colnames
colnames(cov) <- gsub("X282set_","282set_",colnames(cov))
colnames(cov) <- gsub("Goodman.Buckler","",colnames(cov))
colnames(cov) <- gsub("33\\.16","33-16",colnames(cov))
colnames(cov) <- gsub("38\\.11","38-11",colnames(cov))
colnames(cov) <- gsub("A441\\.5","A441-5",colnames(cov))
colnames(cov) <- gsub("CH701\\.30","CH701-30",colnames(cov))
colnames(cov) <- gsub("CI187\\.2","CI187-2",colnames(cov))

# filter exonic ACRs
acr.data$V12[acr.data$V12 == '.'] <- 0
acr.data$V12 <- as.numeric(as.character(acr.data$V12))
acr.data <- subset(acr.data, acr.data$V4 != "exons")
acr.data$V4 <- as.character(acr.data$V4)

# find shared ACRs
rownames(acr.data) <- paste(acr.data$V1,acr.data$V2,acr.data$V3,sep="_")
shared.ids <- intersect(rownames(acr.data), rownames(scatac))
scatac <- scatac[shared.ids,]

# add scATAC + wgs
ids <- as.data.frame(do.call(rbind, strsplit(colnames(scatac),"\\.")))
genotypes <- unique(as.character(ids$V1))
gt <- lapply(genotypes, function(x){
    idx <- which(as.character(ids$V1)==x)
    rowSums(scatac[,idx])
})
gt <- do.call(cbind, gt)
colnames(gt) <- genotypes
row.ids <- intersect(rownames(gt), rownames(wgs))
col.ids <- intersect(colnames(gt), colnames(wgs))
col.ids <- intersect(col.ids, colnames(cov))
gt <- gt[row.ids, col.ids]
wgs <- wgs[row.ids, col.ids]
cov <- cov[row.ids, col.ids]
cov <- as.matrix(cov)
both <- gt+wgs
pav_mat <- both
pav_mat[pav_mat > 0] <- 1
pav_mat <- pav_mat[mixedorder(rownames(pav_mat), decreasing=F),]

# analyze gt ACR accessibility
inbred.data <- inbred.data[col.ids,]
inbred.data$total <- colSums(gt)
inbred.data$open <- colSums(gt > 0)
idf <- subset(inbred.data, inbred.data$open > (0.9*nrow(gt)))
gtf <- gt[,rownames(idf)]
wgs <- wgs[,colnames(gtf)]
cov <- cov[,colnames(gtf)]
both <- both[,colnames(gtf)]

# scale by ACR size
acrs.coords <- as.data.frame(do.call(rbind, strsplit(rownames(wgs),"_")))
acrs.coords$V2 <- as.numeric(as.character(acrs.coords$V2))
acrs.coords$V3 <- as.numeric(as.character(acrs.coords$V3))
acrs.coords$len <- (acrs.coords$V3-acrs.coords$V2)/1000

# scale per kb
gtf <- gtf / acrs.coords$len
wgs <- wgs / acrs.coords$len
both <- both / acrs.coords$len

# identify conserved/non-conserved ACRs
n1 <- rowSums(cov >= 0.9)
n2 <- rowSums(cov <= 0.1)
pdf("ACR_conservation_by_coverage.ri.05.16.23.pdf", width=5, height=5)
plot_cbd(n2, n1, cex=0.8, bwf=50, colP=colorRampPalette(rev(brewer.pal(11, "Spectral"))),
         xlab="N2: n. taxa that align <= 10% of ACR",
         ylab="N1: n. taxa that align >= 90% of ACR", rasterize=T)
dev.off()

# plot phylo scores
phylo <- phylo.data$V2
names(phylo) <- as.character(phylo.data$V1)

# get conserved group
df <- data.frame(n1=n1, n2=n2)
df$phylo <- phylo[rownames(df)]
con <- rownames(df)[df$n1 > 130 & df$n2 < 25]
df$conserved <- ifelse(rownames(df) %in% con, "G1", "G2")
rownames(acr.data) <- paste(acr.data$V1,acr.data$V2,acr.data$V3,sep="_")
type <- acr.data$V4
names(type) <- rownames(acr.data)
df$type <- type[rownames(df)]
rownames(teo) <- paste(teo$V1,teo$V2,teo$V3,sep="_")
shared <- intersect(rownames(df), rownames(teo))
teo <- teo[shared,]

# plot n1/n2 with phylo-scores
df <- df[order(df$phylo, decreasing=F),]
#df <- df[sample(seq(1:nrow(df))),]
pdf("PhyloP_scores.05.16.2023.pdf", width=5, height=5)
plot_raster(df$n2, df$n1, col=colorRampPalette(viridis(10), bias=1.5)(100)[cut(df$phylo, breaks=101)],
     cex=0.8, pch=16, 
     xlab="N2: n. taxa that align <= 10% of ACR",
     ylab="N1: n. taxa that align >= 90% of ACR")
grid(lty=1)
dev.off()

# estimate fold change from genomic dist
prop.group <- prop.table(table(df$conserved, df$type),1)
prop.all <- prop.table(table(df$type))
prop.change <- prop.group %*% diag(1/prop.all)
colnames(prop.change) <- names(prop.all)
lf2 <- log2(prop.change)
pdf("Fold_ACR_change_over_genomic_average.05.16.2023.pdf", width=5, height=5)
pheatmap(lf2, col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         breaks=c(seq(from=min(lf2), to=0, length.out=51), 
                  seq(from=0.00000001, to=max(lf2), length.out=50)))
dev.off()

# plot fold change as bubble
pcc <- melt(log2(prop.change))
pc1 <- melt(table(df$conserved, df$type))
colnames(pcc) <- c("Group","Type","Proportion")
pcc$Total <- pc1$value
pcc$Type <- factor(pcc$Type, levels=c("TSS", "exons", "introns", "promoters", "intergenic"))
pcc$Group <- factor(pcc$Group, levels=c("G1","G2"))
pcc$x <- as.numeric(pcc$Type)
pcc$y <- as.numeric(pcc$Group)
#pcc$area <- rescale((pi*((pcc$Total/nrow(df))^2)), c(2, 5))
pcc$area <- rescale((pi*(abs(pcc$Proportion)^2)), c(2, 5))
cols <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
cols2 <- cols[cut(pcc$Proportion, breaks=c(seq(from=-1.34, to=0, length.out=51),seq(from=0.0000001, to=0.5, length.out=50)))]
pdf("bubble_fold_change_ACRs.pdf", width=5, height=3)
plot(pcc$x, pcc$y, cex=pcc$area, pch=16, col=cols2, ylim=c(0.5,2.5))
grid(lty=1)
dev.off()

# intersect with TEs
te <- read.table("/scratch/apm25309/reference_genomes/Zmays/v5/genomic_features/ACR_TE.ann.bed")
all.te <- read.table(gzfile("/scratch/apm25309/reference_genomes/Zmays/v5/Zm-B73-REFERENCE-NAM-5.0.TE.gff3.gz"))
gw.te <- data.frame(table(all.te$V3))
colnames(gw.te) <- c("TE", "GW_count")
rownames(gw.te) <- gw.te$TE
te$frac <- calcOverlapACR(te)
te <- te[order(te$frac, decreasing=T),]
te <- te[!duplicated(te[,1:3]),]
te <- subset(te, te$frac >= 0.5)
rownames(te) <- paste(te$V1,te$V2,te$V3,sep="_")
te.type <- as.character(te$V7)
names(te.type) <- rownames(te)
df$TE <- te.type[rownames(df)]
df$TE[is.na(df$TE)] <- "nonTE"
df$olTE <- ifelse(df$TE == "nonTE", 0, 1)

# compare proportion ACRs overlap TE between conserved/dispensable
con.te <- table(df$olTE[df$conserved=="G1"])
dis.te <- table(df$olTE[df$conserved=="G2"])
pdf("Pie_chart_TE.05.16.2023.pdf", width=6, height=4)
layout(matrix(c(1:2), nrow=1))
par(mar=c(1,1,1,1))
pie(con.te)
pie(dis.te)
dev.off()

# compare TE classes
c.te <- subset(df, df$olTE==1 & df$conserved=="G1")
n.te <- subset(df, df$olTE==1 & df$conserved=="G2")
a.te <- subset(df, df$olTE==1)
con.te <- table(c.te$TE)
dis.te <- table(n.te$TE)
gw.te$conserved <- con.te[rownames(gw.te)]
gw.te$dispensable <- dis.te[rownames(gw.te)]
gw.te[is.na(gw.te)] <- 0
gw.te$conserved.p <- gw.te$conserved/sum(gw.te$conserved)
gw.te$dispensable.p <- gw.te$dispensable/sum(gw.te$dispensable)
#gw.te$conserved.p <- log2(prop.table(gw.te$conserved)/prop.table(gw.te$GW_count))
#gw.te$dispensable.p <- log2(prop.table(gw.te$dispensable)/prop.table(gw.te$GW_count))
#gw.te$conserved.p[is.infinite(gw.te$conserved.p)] <- 0
#gw.te$dispensable.p[is.infinite(gw.te$dispensable.p)] <- 0


##########################################
##### compare fixed te-ACRs classess #####
##########################################

# add information from the 21 teosinte 
shared.ids <- intersect(rownames(df), rownames(teo))
df <- df[shared.ids,]
teo <- teo[shared.ids,]
df$teo_con <- ifelse((rowSums(teo <= 0.1) = ncol(teo), 0, 1))

# maize fixed te-ACRs
fixed <- subset(df, df$conserved=="G1" & df$olTE==1 & df$teo_con==0)
fixed.props <- data.frame(table(fixed$TE))
f.props <- fixed.props$Freq
names(f.props) <- fixed.props$Var1

# maize te-ACRs
allteACRs.props <- data.frame(table(df$TE[df$TE != "nonTE" & df$conserved=="G1"]))
rownames(allteACRs.props) <- allteACRs.props$Var1
colnames(allteACRs.props) <- c("TE_family", "all.teACRs")
allteACRs.props$maize.fixed.teACRs <- f.props[rownames(allteACRs.props)]
allteACRs.props$maize.fixed.teACRs[is.na(allteACRs.props$maize.fixed.teACRs)] <- 0
allteACRs.props$TE_family <- NULL
teACRs <- as.matrix(allteACRs.props)
teACR.counts <- teACRs
teACRs <- teACRs %*% diag(1/colSums(teACRs))
colnames(teACRs) <- c("all_teACRs", "maize_sp_teACRs")

# test differences
chi.sq <- lapply(seq(1:nrow(teACRs)), function(x){
  spte <- teACR.counts[x,]
  other <- colSums(teACR.counts[!rownames(teACR.counts) %in% rownames(teACR.counts)[x],])
  mm <- rbind(spte, other)
  rownames(mm) <- c(rownames(teACR.counts)[x], "other")
  test <- chisq.test(mm)
  data.frame(TE_family=rownames(teACR.counts)[x], 
             Chisq=as.numeric(test$statistic),
             Pvalue=as.numeric(test$p.value))
})
chi.sq <- do.call(rbind, chi.sq)
chi.sq$log2fc <-  log2((teACRs[,2]+0.01)/(teACRs[,1]+0.01))
chi.sq$signX2 <- chi.sq$Chisq * sign(chi.sq$log2fc)
chi.sq <- chi.sq[order(chi.sq$signX2, decreasing=T),]
chi.sq$FDR <- p.adjust(chi.sq$Pvalue, method="fdr")
rownames(chi.sq) <- chi.sq$TE_family
    
# barplot
pdf("chisquare.pdf", width=5, height=7)
chi.sq <- chi.sq[order(chi.sq$log2fc, decreasing=T),]
cols <- colorRampPalette(brewer.pal(9, "Reds")[2:9])(100)
cols <- cols[cut(chi.sq$log2fc, breaks=seq(from=min(chi.sq$log2fc)-0.01, to=max(chi.sq$log2fc)+0.1, length.out=101))]
par(mar=c(5,15,1,1))
barplot(rev(chi.sq$log2fc), las=2, border=NA, horiz=T, col=rev(cols), names.arg=rev(rownames(chi.sq)))
box()
grid(lty=1)
dev.off()


# plot heatmap of proportions
pdf("teACR_props.pdf", width=5, height=5)
pheatmap(teACRs[rownames(chi.sq),], cluster_row=F, col=viridis(100))
dev.off()




## STARR activity 
df$STARR <- acr.data[rownames(df),]$V12
df$conserved_rate <- df$n1/(df$n1+df$n2)

# fit model
df$Bin_activity1 <- ifelse(df$STARR > 0, 1, 0)
df$Bin_activity2 <- ifelse(df$STARR > 0, df$STARR, NA)
mod <- glm(Bin_activity1~conserved_rate, data=df, family=binomial())

# predict
newx <- seq(from=0,to=1, length.out=500)
preds <- predict(mod, newdata=data.frame(conserved_rate=newx),type="response")
m.val <- mean(df$STARR[df$conserved_rate==1])
preds.new <- rescale(preds, c(0,m.val))

# plot
pdf("conserved_rate_enhancer.pdf", width=5, height=5)
plot_cbd(df$conserved_rate, df$STARR, 
         colP=colorRampPalette(c("grey75",brewer.pal(9, "YlGnBu")), bias=1.5),
         cex=0.8, 
         bwf=25,
         rasterize=T,
         xlab="Conservation rate",
         ylab="Enhancer activity")
lines(newx, preds.new, col="firebrick3", lwd=2)
grid(lty=1)
dev.off()


# MOTIFS MOTIFS MOTIFSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSs
motif <- readRDS("../../step6_co_accessibility/maize_282.v8.3.ALL_CELLs.109k_acrs.motif_matches.rds")
df <- read.table("ACR_Evo_annotation.05.17.2023.txt")
mm <- motif@assays@data$motifMatches
mm.sites <- as.data.frame(motif@rowRanges)
rownames(mm) <- paste(mm.sites$seqnames, mm.sites$start, mm.sites$end, sep="_")
colnames(mm) <- rownames(motif@colData)
mm <- as(mm, "dgCMatrix")

# filter motifs
fm <- mm[rownames(df),]

# run elastic net regression
mod <- cv.glmnet(fm, df$conserved_rate, alpha=0)
l.min <- mod$lambda.min
betas <- coef(mod, s = "lambda.min")
beta <- as.numeric(betas[,1])
names(beta) <- rownames(betas)
beta <- beta[2:length(beta)]
beta <- beta[sort(names(beta))]

# get permuted betas
perms <- mclapply(seq(1:1000), function(z){
    
    # message
    message(" - initialize permutation ", z)
    
    # randomize motif/acr matrix
    rfm <- fm
    colnames(rfm) <- colnames(fm)[sample(ncol(fm))]
    rownames(rfm) <- rownames(rfm)[sample(nrow(fm))]
    rfm <- rfm[rownames(df),]
    rfm <- rfm[,colnames(fm)]
    
    # run model
    rmod <- glmnet(rfm, df$conserved_rate, alpha=0)
    rbetas <- coef(rmod, s = l.min)
    rbeta <- as.numeric(rbetas[,1])
    names(rbeta) <- rownames(rbetas)
    rbeta <- rbeta[2:length(rbeta)]
    rbeta <- rbeta[sort(names(rbeta))]
    return(rbeta)
    
}, mc.cores=24)
perms <- do.call(rbind, perms)
rownames(perms) <- paste0("perms", seq(1:nrow(perms)))
saveRDS(perms, file="permuted_motif_beta.rds")

# compare to observed
perms <- perms[,names(beta)]
p.ave <- colMeans(perms, na.rm=T)
p.sd <- apply(perms, 2, sd, na.rm=T)

# z-score on beta
z.score <- (beta - p.ave)/p.sd
z.score <- z.score[order(z.score, decreasing=T)]
p.value <- pnorm(z.score, lower.tail=F)
q.value <- p.adjust(p.value, method="fdr")

# plot
pdf("Motif_enrichment_Zscore.pdf", width=5, height=5)
cols <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
cols.c <- cols[cut(z.score, breaks=c(seq(from=(min(z.score)), to=0, length.out=51),
                                     seq(from=0.000001, to=max(z.score), length.out=50)))]
plot(z.score, pch=16, cex=1, col=cols.c,
     ylim=c(-1*max(abs(z.score)), max(abs(z.score))),
     xlab="TFBS motif rank",
     ylab="Beta z-score")
grid(lty=1)
dev.off()

# find TFs associated with loss of conservation (not under selection)
z.score <- z.score[order(z.score, decreasing=F)]
lp.value <- pnorm(z.score, lower.tail=T)
lq.value <- p.adjust(lp.value, method="fdr")


# test motif enhancer activity
emod <- cv.glmnet(fm, df$phylo, alpha=0)
ebetas <- coef(emod, s = "lambda.min")
ebeta <- as.numeric(ebetas[,1])
names(ebeta) <- rownames(ebetas)
ebeta <- ebeta[2:length(ebeta)]
ebeta <- ebeta[sort(names(ebeta))]


##########################################################################################################################################
# normalized matrices
gtf.norm <- cpm(gtf, log=F)
wgs.norm <- cpm(wgs, log=F)
both.norm <- cpm(both, log=F)
cov.pav <- cov
cov.pav[cov.pav < 0.5] <- 0
cov.pav[cov.pav > 0] <- 1

# count the number of sites with no WGS and + ATAC data
dispen <- lapply(colnames(both.norm), function(z){
    ifelse(both.norm[,z] < 1 & cov[,z] < 0.5, 0, 1)
})
dispen <- do.call(cbind, dispen)
colnames(dispen) <- colnames(both.norm)
acr.cnt <- rowSums(dispen==0)
ind.cnt <- colSums(dispen==0)
head(ind.cnt)
table(acr.cnt)

# plot
pdf("PAV_ACRs.barplot_ACRs.dispensable.pdf", width=6, height=6)
hist(acr.cnt, col="grey", breaks=2*max(acr.cnt), border=NA)
dev.off()

# compare dispensable/non-dispensable between tropical and non-tropical lines
inbred.data <- inbred.data[colnames(gtf),]
trop <- dispen[,rownames(inbred.data)[inbred.data$type=="tropical"]]
nontrop <- dispen[,rownames(inbred.data)[inbred.data$type=="nontropical"]]
trop.frac <- rowMeans(trop)
nontrop.frac <- rowMeans(nontrop)
dif <- nontrop.frac - trop.frac

# find accessible/inaccessible threshold
bgt <- gtf
use.zTPM <- T
if(use.zTPM){
    thresh <- -2
    zgtf <- as.matrix(zFPKM(as.data.frame(log2(gtf+1))))
    bgt[zgtf > thresh] <- 1
    bgt[zgtf <= thresh] <- 0
}else{
    thresh <- 2
    bgt[bgt < thresh] <- 0
    bgt[bgt > 0] <- 1
}
freq <- rowSums(bgt)

trop <- bgt[,rownames(inbred.data)[inbred.data$type=="tropical"]]
nontrop <- bgt[,rownames(inbred.data)[inbred.data$type=="nontropical"]]
trop.frac <- rowMeans(trop)
nontrop.frac <- rowMeans(nontrop)
dif <- nontrop.frac - trop.frac

# plot
pdf("ACR_open_closed_frequency.pdf", width=6, height=6)
plot(density(freq, from=0, to=171, bw=1), main="Accessibility Frequency", xlab="# of accessible genotypes")
dev.off()

# save data
bed <- read.table("/scratch/apm25309/reference_genomes/Zmays/v5/genomic_features/celltype_genotype.ACRs.109k.features.bed")
ff <- as.data.frame(freq)
rownames(bed) <- paste(bed$V1,bed$V2,bed$V3,sep="_")
ff$ACRtype <- ifelse(ff$freq >= ceiling(ncol(gtf)*0.95), "Core", "Variable")
bed$ACRtype <- ff[rownames(bed),]$ACRtype
bed$num_geno_open <- ff[rownames(bed),]$freq
gtf <- gtf[rownames(bed),]
bgt <- bgt[rownames(bed),]
bed$ave_Tn5_density <- unlist(lapply(seq(1:nrow(gtf)), function(z){
    nums <- gtf[z,]
    nums <- nums[as.logical(bgt[z,])]
    mean(nums)
}))
props <- prop.table(table(bed$ACRtype, bed$V4), 1)
cc <- props[1,]
vv <- props[2,]

# plot pie chart
pdf("Pie_dispensable_ACRs_genomic_dist.pdf", width=10, height=5)
layout(matrix(c(1:2), nrow=1, ncol=2, byrow=T))
pie(cc)
pie(vv)
dev.off()

# id PAV
num.missing <- colSums(both==0)
num.missing <- num.missing[order(num.missing, decreasing=T)]
shared <- intersect(names(num.missing), rownames(inbred.data))
inbred.data <- inbred.data[shared,]
ii <- as.matrix(inbred.data[,c(2:6)])
pdf("Subpopulation_sorted_ACR_dispensable.pdf", width=12, height=4)
image(t(ii), col=colorRampPalette(c("white", "firebrick4"))(100))
dev.off()

pdf("PAV_ACRs.barplot_genotype.pdf", width=12, height=6)
barplot(log10(num.missing+1),las=2, border=NA)
dev.off()

# count number of genotypes an ACR has disappeared in
g.counts <- rowSums(both==0)
pdf("PAV_ACRs.barplot_ACRs.pdf", width=6, height=6)
hist(g.counts, col="grey", border=NA, breaks=60)
dev.off()


# make df
df <- as.data.frame(do.call(rbind, strsplit(names(g.counts),"_")))
df$V4 <- g.counts
df$V2 <- as.numeric(as.character(df$V2))
df$V3 <- as.numeric(as.character(df$V3))
df$V1 <- as.character(df$V1)
df <- df[order(df$V1,df$V2,decreasing=F),]
bed <- read.table("/scratch/apm25309/reference_genomes/Zmays/v5/genomic_features/celltype_genotype.ACRs.109k.features.bed")
bed$V5 <- df$V4
core <- subset(bed, bed$V5 <= 1)
disp <- subset(bed, bed$V5 > 1)
p.core <- prop.table(table(core$V4))
p.disp <- prop.table(table(disp$V4))

# plot pie chart
pdf("Pie_dispensable_ACRs_genomic_dist.pdf", width=10, height=5)
layout(matrix(c(1:2), nrow=1, ncol=2, byrow=T))
pie(p.core)
pie(p.disp)
dev.off()

# subgenome data
sg <- read.table("maize_subgenome_homeologs.ACRs.bed")
sg <- sg[sg$V12 <= 2000,]
sg <- sg[order(sg$V12, decreasing=F),]
sg$acrID <- paste(sg$V8,sg$V9,sg$V10,sep="_")
core$acrID <- paste(core$V1,core$V2,core$V3,sep="_")
disp$acrID <- paste(disp$V1,disp$V2,disp$V3,sep="_")
sg <- sg[!duplicated(sg$acrID),]
subg <- sg$V7
names(subg) <- sg$acrID
core$subgenome <- subg[core$acrID]
disp$subgenome <- subg[disp$acrID]
core.sub <- prop.table(table(core$subgenome))
disp.sub <- prop.table(table(disp$subgenome))

pdf("Pie_dispensable_ACRs_subgenome.pdf", width=10, height=5)
layout(matrix(c(1:2), nrow=1, ncol=2, byrow=T))
pie(core.sub)
pie(disp.sub)
dev.off()

# TE data
te <- read.table("/scratch/apm25309/reference_genomes/Zmays/v5/genomic_features/test.1")
te <- te[!duplicated(te[,1:3]),]
rownames(te) <- paste(te$V1,te$V2,te$V3,sep="_")
rownames(core) <- paste(core$V1,core$V2,core$V3,sep="_")
rownames(disp) <- paste(disp$V1,disp$V2,disp$V3,sep="_")
te.id <- as.character(te$V11)
names(te.id) <- rownames(te)
core$TE <- te.id[rownames(core)]
disp$TE <- te.id[rownames(disp)]
core$TE[is.na(core$TE)] <- 0
disp$TE[is.na(disp$TE)] <- 0
core$isTE <- ifelse(core$TE==0, 0, 1)
disp$isTE <- ifelse(disp$TE==0, 0, 1)

# intergenic only
i.core <- subset(core, core$V4=="intergenic")
i.disp <- subset(disp, disp$V4=="intergenic")

pdf("Pie_dispensable_ACRs_TE_Prop.pdf", width=10, height=5)
layout(matrix(c(1:2), nrow=1, ncol=2, byrow=T))
pie(prop.table(table(i.core$isTE)))
pie(prop.table(table(i.disp$isTE)))
dev.off()

# prop type TE
i.core.TE <- prop.table(table(i.core$TE[i.core$TE != 0]))*100
i.disp.TE <- prop.table(table(i.disp$TE[i.disp$TE != 0]))*100
ids <- c(names(i.core.TE), names(i.disp.TE))
ids.u <- unique(ids)
dff <- data.frame(coreTE=i.core.TE[ids.u], dispTE=i.disp.TE[ids.u], row.names=ids.u)
dff$coreTE.Var1 <- NULL
dff$dispTE.Var1 <- NULL
dff <- as.matrix(dff)
dff[is.na(dff)] <- 0
colnames(dff) <- c("core", "dispensable")

pdf("TE_proportions.intergenicTEs.pdf", width=6, height=7)
pheatmap(dff, col=viridis(100))
dev.off()
pheatmap(dff, col=viridis(100))
