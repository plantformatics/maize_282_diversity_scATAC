## multi-acr caQTL v2 analysis ##

# libraries
library(vioplot)

# functions
locateVar <- function(x){
  
  snp <- as.data.frame(do.call(rbind, strsplit(x$variant_id,"_")))
  snp$V2 <- as.numeric(snp$V2)
  acr1 <- as.data.frame(do.call(rbind, strsplit(x$acr1,"_")))
  acr2 <- as.data.frame(do.call(rbind, strsplit(x$acr2,"_")))
  acr1$V2 <- as.numeric(acr1$V2)
  acr1$V3 <- as.numeric(acr1$V3)
  acr2$V2 <- as.numeric(acr2$V2)
  acr2$V3 <- as.numeric(acr2$V3)
  acr1$inACR <- ifelse(snp$V2 >= acr1$V2 & snp$V2 <= acr1$V3, 1, 0)
  acr2$inACR <- ifelse(snp$V2 >= acr2$V2 & snp$V2 <= acr2$V3, 1, 0)
  x$varIn <- ifelse(acr1$inACR==1, 1, 2)
  return(x)
  
}

# load data
a <- readRDS("long_range_caQTL_interactions.rds")
b <- read.table("/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/step10_caQTL_celltypes/celltype_genotype.ACRs.109k.features.STARR.bed")
caqtl <- read.table("../../step10_caQTL_celltypes/fastQTL_testing/celltype_results_09_2023/All_caQTL_ACRs.txt")

# pinpoint variant
a <- locateVar(a)

# add info
rownames(b) <- paste(b$V1,b$V2,b$V3,sep="_")
b$V12[b$V12=='.'] <- 0
b$V12 <- as.numeric(b$V12)
starr <- b$V12
names(starr) <- rownames(b)
type <- b$V4
names(type) <- rownames(b)
a$acr1STARR <- starr[a$acr1]
a$acr2STARR <- starr[a$acr2]
a$acr1Type <- type[a$acr1]
a$acr2Type <- type[a$acr2]

focal <- data.frame(starr=ifelse(a$varIn==1, a$acr1STARR, a$acr2STARR),
                    type=ifelse(a$varIn==1, a$acr1Type, a$acr2Type))

nonfocal <- data.frame(starr=ifelse(a$varIn==2, a$acr1STARR, a$acr2STARR),
                       type=ifelse(a$varIn==2, a$acr1Type, a$acr2Type))
qtl <- unique(caqtl$V1)
                    
                    
focal.props <- table(focal$type)
back.props <- table(b$V4[b$V11=="tested"])
qtl.props <- prop.table(table(b$V4[rownames(b) %in% qtl]))
mat <- matrix(as.numeric(c(focal.props["intergenic"], nrow(focal), back.props["intergenic"], length(b$V11[b$V11=="tested"]))), nrow=2, byrow=T)
chisq.test(mat)$p.value
rownames(both) <- c("lrACRs", "allACRs")

pdf("distribution_lrACR_type.pdf", width=5, height=4)
barplot(both, beside=T)
grid(lty=1)
box()
dev.off()

# compare enhancer activity
nonlr <- data.frame(starr=starr[qtl],
                    type=type[qtl])


distal.focal <- subset(focal, focal$type=="intergenic")
distal.nonlr <- subset(nonlr, nonlr$type=="intergenic")
distal.nonqtl <- subset(b, b$V4=="intergenic" & b$V11!="tested")

pdf("starr_activity.pdf", width=5, height=5)
boxplot(distal.focal$starr, distal.nonlr$starr, distal.nonqtl$V12, outline=F)
dev.off()