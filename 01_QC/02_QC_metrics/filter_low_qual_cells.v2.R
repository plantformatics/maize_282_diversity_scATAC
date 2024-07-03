# remove bad cells #

# load args
args <- commandArgs(T)
if(length(args) != 2){stop("Rscript filter_low_qual_cells.R <input.rds> <preix>")}
input.dat <- as.character(args[1])
prefix <- as.character(args[2])

# load libraries
library(Socrates)
library(qlcMatrix)

# load data
message(" - loading data")
a <- readRDS(input.dat)
prefilt <- read.table("All_pools.clusters.tsv")

# align ids with counds
shared <- intersect(rownames(a$meta), colnames(a$counts))
a$counts <- a$counts[,shared]
a$meta <- a$meta[shared,]

# generate stats
a$meta <- a$meta[order(a$meta$nSites, decreasing=T),]
a$meta$pTSS <- a$meta$tss/a$meta$total
a$meta$FRiP <- a$meta$acrs/a$meta$total
a$meta$pOrg <- a$meta$ptmt/a$meta$total
a$meta$qc_check <- ifelse(rownames(a$meta) %in% rownames(prefilt), 1, 0)

# set initial thresholds
message(" - setting filters")
a$meta$tss_z <- (a$meta$pTSS - mean(a$meta$pTSS[a$meta$qc_check==1]))/sd(a$meta$pTSS[a$meta$qc_check==1])
a$meta$acr_z <- (a$meta$FRiP - mean(a$meta$FRiP[a$meta$qc_check==1]))/sd(a$meta$FRiP[a$meta$qc_check==1])
a$meta$sites_z <- (log10(a$meta$nSites) - mean(log10(a$meta$nSites[a$meta$qc_check==1])))/sd(log10(a$meta$nSites[a$meta$qc_check==1]))
a$meta$tss_z[is.na(a$meta$tss_z)] <- -10
a$meta$acr_z[is.na(a$meta$acr_z)] <- -10
a$meta$sites_z[is.na(a$meta$sites_z)] <- -10 
a$meta$qc_check <- ifelse(a$meta$tss_z < -2 | a$meta$pTSS < 0.2, 0, 
                          ifelse(a$meta$acr_z < -2 | a$meta$FRiP < 0.2, 0,
				ifelse(a$meta$sites_z < -2, 0, a$meta$qc_check)))

# parse
message(" - parsing initial boundaries")
good.cells <- rownames(subset(a$meta, a$meta$qc_check==1))
bad.cells <- rownames(subset(a$meta, a$meta$qc_check==0))
gg <- a$counts[,colnames(a$counts) %in% good.cells]
gg <- gg[,Matrix::colSums(gg > 0) > 100]
a$meta$qc_check <- ifelse(rownames(a$meta) %in% colnames(gg), 1, 0)
bb <- a$counts[,! colnames(a$counts) %in% colnames(gg)]
sites <- Matrix::rowMeans(gg > 0)
sites <- sites[order(sites, decreasing=T)]
num.sites <- min(c(max(a$meta$nSites), 250000))
if(length(sites) < num.sites){
	num.sites <- length(sites)
}
gg <- gg[names(sites)[1:num.sites],]
gg <- gg[,Matrix::colSums(gg) > 0]
bb <- bb[rownames(gg),]
bb <- bb[,Matrix::colSums(bb) > 0]
shared <- intersect(rownames(gg), rownames(bb))
gg <- gg[shared,]
bb <- bb[shared,]
gg <- gg[,Matrix::colSums(gg) > 0]
bb <- bb[,Matrix::colSums(bb) > 0]

# Do not use more than 250,000 'bad' cells
if(ncol(bb) > 250000){
    top <- Matrix::colSums(bb)
    top <- top[order(top, decreasing=T)]
    top <- names(top)[5001:250000] # skip cells at the boundary
}else{
    num.cells <- ncol(bb)
    top <- Matrix::colSums(bb)
    top <- top[order(top, decreasing=T)]
    if(num.cells > 100000){
        top <- names(top)[5001:100000]
    }else{
        top <- names(top)
    }
}

# clean up ref
bb <- bb[,top]
bb <- bb[Matrix::rowSums(bb)>0,]
shared.sites <- intersect(rownames(bb), rownames(gg))
bb <- bb[shared.sites,]
gg <- gg[shared.sites,]
bb <- bb[,Matrix::colSums(bb)>0]
gg <- gg[,Matrix::colSums(gg)>0]

# make references
num.good <- 1000
if(ncol(gg) < num.good){
	num.good <- ncol(gg)
}
top.gg <- Matrix::colSums(gg)
top.gg <- top.gg[order(top.gg, decreasing=T)]
top.gg.ids <- names(top.gg)[1:num.good]

# function
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

message(" - normalizing distributions and creating references")
sub.counts <- a$counts[,c(colnames(gg), colnames(bb))]
sub.counts <- sub.counts[rownames(gg),]
all.res <- tfidf(list(counts=sub.counts), doL2=T)$residuals
bb.norm <- all.res[,colnames(bb)]
gg.norm <- all.res[,colnames(gg)]

# pick sites
res.ave <- Matrix::rowMeans(gg.norm)
res.res <- RowVar(gg.norm)
resis <- loess(res.res~res.ave)$residuals
names(resis) <- rownames(gg.norm)
#resis <- resis[order(resis, decreasing=T)]
top.sites <- names(resis[resis > 0])
bb.norm <- bb.norm[top.sites,]
gg.norm <- gg.norm[top.sites,]

# make references
bad.ref <- Matrix::rowMeans(bb.norm)
good.ref <- Matrix::rowMeans(gg.norm) #gg.norm[,top.gg.ids]
bad.ref <- Matrix(matrix(c(bad.ref / (sqrt(sum(bad.ref^2)))), ncol=1), sparse=T)
good.ref <- Matrix(matrix(c(good.ref / (sqrt(sum(good.ref^2)))),ncol=1), sparse=T)


# check each cell against ref
message(" - estimating correlations")
b.ref <- corSparse(gg.norm, bad.ref)
g.ref <- corSparse(gg.norm, good.ref)
b.bad <- corSparse(bb.norm, bad.ref)
g.bad <- corSparse(bb.norm, good.ref)
refs <- data.frame(cbind(b.ref, g.ref))
bads <- data.frame(cbind(b.bad, g.bad))
colnames(refs) <-c("bREF", "gREF")
colnames(bads) <-c("bREF", "gREF")
refs$call <- ifelse(refs$bREF > refs$gREF, 0, 1)
bads$call <- ifelse(bads$bREF > bads$gREF, 0, 1)
rownames(refs) <- colnames(gg.norm)
rownames(bads) <- colnames(bb.norm)
all.refs <- rbind(refs,bads)
shared <- intersect(rownames(a$meta), rownames(all.refs))
meta <- a$meta[shared,]
all.refs <- all.refs[shared,]
all.refs <- cbind(meta, all.refs)
nonrefs <- a$meta[!rownames(a$meta) %in% rownames(all.refs),]
nonrefs$bREF <- NA
nonrefs$gREF <- NA
nonrefs$call <- NA
test <- rbind(all.refs, nonrefs)

# get new call
test$dif <- test$gREF-test$bREF
top <- test$dif[test$qc_check==1]
bottom <- test$dif[test$qc_check==0]
cut.off <- median(bottom, na.rm=T)
test$call <- ifelse(test$dif > cut.off & test$qc_check == 1, 1, 0)
test$dif <- NULL
write.table(test, file=paste0(prefix,".updated_metadata.v3.txt"), quote=F, row.names=T, col.names=T, sep="\t")
