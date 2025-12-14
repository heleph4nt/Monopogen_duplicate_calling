library(e1071)        # Provides SVM implementation used for classification
library(reshape2)    # Used for data reshaping (long ↔ wide formats)

# Source helper functions for SVM-based phasing
source("/rsrch3/scratch/bcb/jdou1/scAncestry/Monopogen/src/svm_phasing_sub.R")

# Read command-line arguments
args = commandArgs(trailingOnly=TRUE)

chr <- args[1]                      # Chromosome identifier
loaddata <- args[2]                 # Flag: whether to load raw data or an RDS image
dis_limit_physical <- as.numeric(args[3])  # Physical distance threshold (read length)
dis_limit_genetic <- as.numeric(args[4])   # Genetic distance threshold (LD-based)



print(chr)

# Case 1: load raw data and prepare inputs
if(loaddata==1){
	print("Loading data now...")

	# Load single-cell genotype matrix (cells × SNVs)
	mat <- read.table(paste0(file="../germline/chr",chr,".gl.filter.hc.cell.mat.gz"))

	# Load per-SNV feature information (used for SVM QC)
	feature_info <- read.csv(paste0(file="../debug/chr",chr,".gl.feature.info"),header=F)

	# Create a unique SNV identifier (chr:pos)
	mat$V2 <- paste0(mat$V1,":",mat$V2)

	# Extract meta-information columns
	meta <- mat[,c("V1","V2","V3","V4")]
	meta <- as.data.frame(meta)

	# Initialize flags for variant origin (currently unused)
	meta$in_normal <- 0
	meta$in_tumor <- 0
	meta$somatic <- 0

	# Compute B-allele frequency from genotype matrix
	meta$baf <- frq(mat)

	# Select SVM-relevant features for SNVs present in meta
	svm_info <- feature_info[feature_info$V1%in%meta$V2,]
	svm_dt <- getFeatureInfo(svm_info=svm_info)

	# Keep only SNVs present in both genotype matrix and feature table
	overlap <- intersect(meta$V2, svm_dt$id)

	meta_qc <- meta[meta$V2%in%overlap,]
	mat_qc <- mat[meta$V2%in%overlap,]
	svm_dt_qc <- svm_dt[svm_dt$id%in%overlap,]

	# Add BAF as an SVM feature
	svm_dt_qc$BAF <- meta_qc$baf

	# Save the full workspace for reuse
	save.image(file=paste0(chr,".input.image.RDS"))
}

# Case 2: reload previously saved workspace
if(loaddata==2){
	load(file=paste0(chr,".input.image.RDS"))
}

# Re-source SVM utilities (safety)
source("/rsrch3/scratch/bcb/jdou1/scAncestry/Monopogen/src/svm_phasing_sub.R")
print(args)

print("start to run SVM and haplotype filtering")
print("new version ")

# Identify blocks of nearby SNVs likely belonging to the same haplotype
mutation_block <- SNV_block(summary=meta_qc)

# Prepare SVM input labels
svm_in <- SVM_prepare(mutation_block)

# Train SVM classifier to filter low-quality / spurious SNVs
svm_out <- SVM_train(label =svm_in, data=svm_dt_qc, downsampling=2)

# Define SVM-based filtering thresholds
p_lower <- sort(svm_out$test$POS)[floor(nrow(svm_out$test)*0.1)]
p_upper <- sort(svm_out$test$POS)[floor(nrow(svm_out$test)*0.95)]
filter1 <- (svm_out$test$POS>p_lower)
filter2 <- (svm_out$test$baf>0.1 )

# Retain high-confidence SNVs for phasing
svm_pass <- svm_out$test[filter1 & filter2,]

# Perform haplotype phasing using distance constraints
svm_phasing <- Phasing(
	x=mat_qc,
	snv_meta=svm_pass,
	dis_limit=dis_limit_genetic,
	readlength=dis_limit_physical
)

# Store results
svm_out$phasing <- svm_phasing
svm_out$feature <- svm_dt_qc

# Save final SVM + phasing output
saveRDS(svm_out,file=paste0("svm",chr,".out.RDS"))



# Explicit variable assignments (debug / interactive use)
x=mat_qc
snv_meta=svm_pass
dis_limit=dis_limit_genetic
readlength=dis_limit_physical
