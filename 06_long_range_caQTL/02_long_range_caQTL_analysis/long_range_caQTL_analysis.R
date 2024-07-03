## compare coACR QTL and caQTL ##

# libraries
library(gtools)

# functions
varPos <- function(x){
  
  acr <- as.data.frame(do.call(rbind, strsplit(x$phenotype_id, "_")))
  snp <- as.data.frame(do.call(rbind, strsplit(x$variant_id, "_")))
  acr$V2 <- as.numeric(acr$V2)
  acr$V3 <- as.numeric(acr$V3)
  snp$V2 <- as.numeric(snp$V2)
  x$inACR <- ifelse(snp$V2 >= acr$V2 & snp$V2 <= acr$V3, 1, 0)
  return(x)
  
}
reformatInts <- function(x){
  
  # split acr 1
  acr1 <- as.data.frame(do.call(rbind, strsplit(x$acr1, "_")))
  acr2 <- as.data.frame(do.call(rbind, strsplit(x$acr2, "_")))
  acr1$V2 <- as.numeric(as.character(acr1$V2))
  acr1$V3 <- as.numeric(as.character(acr1$V3))
  acr2$V2 <- as.numeric(as.character(acr2$V2))
  acr2$V3 <- as.numeric(as.character(acr2$V3))
  ave1 <- as.integer((acr1$V2+acr1$V3)/2)
  ave2 <- as.integer((acr2$V2+acr2$V3)/2)
  acrA <- ifelse(ave1 < ave2, x$acr1, x$acr2)
  acrB <- ifelse(ave1 > ave2, x$acr1, x$acr2)
  return(data.frame(acr1=acrA, acr2=acrB, variant_id=x$variant_id))
}

# single locus QTL
a <- read.table("Merged_Genotype.tensorQTL_perms.cis_qtl.nominal_pass.txt", header=T)
a <- varPos(a)

# co-access data
co <- read.table("/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/step6_co_accessibility/coaccessible_ACR_metadata.05.14.2024.txt", header=T)

# iterate over each cell type
outs <- lapply(unique(a$celltype), function(z){
  message(" - identifying genetic interactions of cell type ",z)
  ac <- subset(a, a$celltype==z)
  cnts <- table(ac$variant_id)
  aa <- ac[ac$variant_id %in% names(cnts[cnts>1]),]
  aa$associated <- 1
  aa <- aa[mixedorder(aa$variant_id),]
  keep.snps <- subset(aa, aa$inACR==1)
  var.keep <- unique(keep.snps$variant_id)
  aa <- aa[aa$variant_id %in% var.keep,]
  acrc <- as.data.frame(do.call(rbind, strsplit(aa$phenotype_id, "_")))
  aa <- cbind(aa, acrc)
  
  # iterate over each caQTL-ACR pair
  its <- 0
  ints <- lapply(unique(aa$variant_id), function(x){
    its <<- its + 1
    if((its %% 1000)==0){message(" - iterated over ", its, " caQTL ...")}
    df <- subset(aa, aa$variant_id==x)
    facr <- unique(df$phenotype_id[as.logical(df$inACR)])
    nacr <- unique(df$phenotype_id[!as.logical(df$inACR)])
    combs <- expand.grid(facr, nacr)
    colnames(combs) <- c("acr1", "acr2")
    combs$acr1 <- as.character(combs$acr1)
    combs$acr2 <- as.character(combs$acr2)
    combs$variant_id <- x
    return(combs)
  })
  ints <- do.call(rbind, ints)
  rf <- reformatInts(ints)
  rf <- rf[mixedorder(rf$acr1),]
  rf$phenotype_id <- paste0(rf$acr1,"-",rf$acr2)  
  rf$celltype <- z
  return(rf)
})
rf <- do.call(rbind, outs)

# find shared co-accessible ACRs / long-range caQTL
shared.co <- intersect(rf$phenotype_id, co$phenotype_id)

# fine-mapped
b <- read.table("../../step10_caQTL_celltypes/Finemapping/fine_mapped_SNVs.cs95.100bp_ACR.txt", header=T)
shared <- intersect(rf$variant_id, b$snvID)
fm.rf <- rf[rf$variant_id %in% shared,]

# hic
loops <- read.table("../coACR_supported_hic_hichip.txt", header=T)
loops$phenotype_id <- paste0(loops$chr1,"_",loops$s1,"_",loops$e1,"-",loops$chr2,"_",loops$s2,"_",loops$e2)
shared.l <- intersect(loops$phenotype_id, rf$phenotype_id)
shared.fm <- intersect(loops$phenotype_id, fm.rf$phenotype_id)

# coACR QTL
b <- read.table("Merged_Genotype.coACR_tensorQTL_perms.cis_qtl.01.17.2024.txt", header=T)
ids <- read.table("ACR_id_conv.txt")

# compare caQTL-pair with coACR caQTL
shared <- intersect(rf$phenotype_id, b$phenotype_id)
not.shared <- unique(rf$phenotype_id[!rf$phenotype_id %in% shared])
b.shared <- b[b$phenotype_id %in% shared,]
b.notshared <- b[!b$phenotype_id %in% shared,]