## process 282 maize leaf data ##

# load libraries
library(Socrates)

# arguments
args <- commandArgs(T)
if(length(args) != 4){stop("Rscript QC_scATAC_data.R <bed> <prefix> <gff/gtf> <chr.fai>")}

# load data
bed <- as.character(args[1])
out <- as.character(args[2])
ann <- as.character(args[3])
chr <- as.character(args[4])

# load objects
obj <- loadBEDandGenomeData(bed, ann, chr, attribute="Parent")

# count organellar reads
obj <- countRemoveOrganelle(obj, org_scaffolds=c("scaf_34", "scaf_45", "scaf_23", "scaf_66", "scaf_106", "scaf_113", "scaf_36", "Pt", "Mt"), remove_reads=T)

# call ACRs
obj <- callACRs(obj, genomesize=1.6e9, 
                shift= -75, 
                extsize=150,
                fdr=0.1,
                output=paste0(out,"_peaks"), 
                tempdir=paste0(out, '_macs2_temp'), 
                verbose=T)

# build metadata
obj <- buildMetaData(obj,
                     tss.window=2000, 
                     verbose=TRUE)

# generate sparse matrix
obj <- generateMatrix(obj,
                      filtered=F, 
                      windows=500, 
                      peaks=F, 
                      verbose=T)

# convert to Socrates format for downstream analysis. 
soc.obj <- convertSparseData(obj, verbose=T)

# save QC object
saveRDS(soc.obj, file=paste0(out,".raw.soc.rds"))

