#!/usr/bin/env Rscript

# Original Author: Prof. Folkers
# Modified for Snakemake integration

# Command line arguments
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript ballgown.R <workdir> <metadata_file> <sample_paths_file>", call.=FALSE)
}

workdir <- args[1]           # This will be dge_dir for outputs
mdfile <- args[2]            # Path to metadata file
sample_paths_file <- args[3] # Path to the file containing list of ballgown sample directories

# Set working directory (where outputs will be saved)
setwd(workdir)
getwd() # Print current working directory for verification

# Load libraries
library(ballgown)
library(genefilter)
# library(dplyr) # Not used in original script, but useful for data manipulation
# library(ggplot2) # For plotting, if added later

# Read metadata
print(paste("Reading metadata from:", mdfile))
metadata <- read.csv(mdfile)

# Read sample paths from the provided file
print(paste("Reading sample paths from:", sample_paths_file))
sample_vector <- scan(sample_paths_file, what="", sep="\n")
print("Sample paths to be loaded by Ballgown:")
print(sample_vector)

# Check if sample paths exist - Ballgown will error if not, but an early check can be useful
# for (sample_path in sample_vector) {
#   if (!file.exists(sample_path)) {
#     stop(paste("Error: Sample path not found:", sample_path), call.=FALSE)
#   }
# }

# Load ballgown data structure
# sampletype = "condition" (this should match a column in your metadata file)
# timecourse = FALSE (set to TRUE if you have time course data)
# Note: Ballgown expects the 'ids' column in metadata to match the directory names of the samples.
# The sample_vector should provide paths like ".../ballgown_inputs/sampleA", ".../ballgown_inputs/sampleB", etc.
# Ballgown extracts the sample names (sampleA, sampleB) from these paths.
print("Loading Ballgown data structure...")
bg <- ballgown(samples = sample_vector, dataDir = dirname(sample_vector[1]), meas='all', pData=metadata)
# dataDir is the parent directory of sample folders. dirname(sample_vector[1]) assumes all sample folders share the same parent.

# Filter low-abundance genes
# Here, "low-abundance" means an average of <1 read per base in at least one condition.
# This is a common filtering step to improve statistical power.
print("Filtering low-abundance genes/transcripts...")
bg_filt <- subset(bg, "rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Identify transcripts that are differentially expressed
# sampletype should match a column name in your phenodata (metadata) file.
# E.g., if your metadata has a column named "condition" distinguishing groups.
# Ensure the column name used here exists in your metadata.csv.
# The original script had "group" hardcoded. Let's try to infer it or make it a parameter.
# For now, assuming the first column of metadata after "ids" is the condition/group.
# This is a common convention but might need adjustment.
# A better way: require metadata to have a specific column name e.g., "condition".
# Let's assume metadata has a column named "condition" for group comparison.
if (!"condition" %in% colnames(pData(bg_filt))) {
    stop("Error: Metadata file must contain a column named 'condition' for DGE analysis.", call.=FALSE)
}
print("Performing DGE analysis on transcripts...")
results_transcripts <- stattest(bg_filt, feature='transcript', meas='FPKM', covariate='condition', getFC=TRUE)

# Identify exons that are differentially expressed
# print("Performing DGE analysis on exons...") # This was in original, but might not be primary output
# results_exons <- stattest(bg_filt, feature='exon', meas='FPKM', covariate='condition', getFC=TRUE)


# Add gene names to results
# This assumes your GTF had gene_name attributes. Ballgown TRIES to extract these.
# If not, gene names might be NAs.
print("Adding gene names to transcript results...")
results_transcripts <- data.frame(geneNames=ballgown::geneNames(bg_filt), results_transcripts)


# Sort results by q-value
results_transcripts <- results_transcripts[order(results_transcripts$qval),]

# Write results to CSV files
print(paste("Writing transcript DGE results to:", paste(workdir, "/dge_transcripts.csv", sep="")))
write.csv(results_transcripts, paste(workdir,"/dge_transcripts.csv",sep=""), row.names=FALSE)

# print(paste("Writing exon DGE results to:", paste(workdir, "/dge_exons.csv", sep="")))
# write.csv(results_exons, paste(workdir,"/dge_exons.csv",sep=""),row.names=FALSE) # If exon results are needed

# Extract and write transcript expression values
trans_expr <- texpr(bg_filt, 'all')
print(paste("Writing transcript expression matrix to:", paste(workdir, "/trans_expr.csv", sep="")))
write.csv(trans_expr, paste(workdir,"/trans_expr.csv",sep=""), row.names=TRUE)


# Example of how to plot transcripts (commented out for non-interactive Snakemake runs)
# gene_of_interest <- "MSTRG.105" # Replace with an actual gene ID from your data
# if (gene_of_interest %in% results_transcripts$id) {
#   print(paste("Plotting transcript:", gene_of_interest))
#   plotTranscripts(gene=gene_of_interest, gown=bg_filt, samples=sampleNames(bg_filt),
#                   meas='FPKM', colorby='condition',
#                   main=paste('transcripts from gene', gene_of_interest))
# }


# Example of plotting average transcript expression (commented out)
# plotMeans('MSTRG.105', bg_filt, groupvar='condition', meas='FPKM', colorby='condition')

print("Ballgown DGE script finished.")
