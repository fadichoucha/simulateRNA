#!/usr/bin/env Rscript

suppressWarnings(library(dplyr))
suppressWarnings(library(polyester))
suppressWarnings(library(Biostrings))

#-------------------------------------------------------------
gene_sets_file <- '../output/gene_sets_with_indices.tsv'
output_logfc_table <- '../output/logfc_table.tsv'
fasta_file <- '../data/genome/gencode.v38.transcripts.fa'
output_dir <- '../output/simulated_rnaseq_data'
set_seed <- 123
STN <- 0.1
num_reps <- 5
reads_per_transcript_count <- 300

# Define fold change values for each gene set
fc_values <- list(
  M47362 = 0.5,  
  M47501 = 0.5,    
  M47365 = 1,     
  M47478 = 1,  
  M47412 = 3,     
  M47419 = 5,     
  M47422 = 5,     
  M47425 = 10,    
  M47450 = 10,  
  M47452 = 20   
)

#-------------------------------------------------------------

# Read the TSV file
gene_sets <- read.table(gene_sets_file, header=TRUE, sep='\t')

# Add LogFC value to geneset
set.seed(set_seed)

logfc_table <- gene_sets %>%
  group_by(gene_set) %>%
  mutate(fc = runif(n(), fc_values[[unique(gene_set)]] - STN, fc_values[[unique(gene_set)]] + STN))

write.table(logfc_table, file = output_logfc_table, sep = '\t', row.names = FALSE, quote = FALSE)

# Load the FASTA sequences
fasta_sequences <- readDNAStringSet(fasta_file)

# Create a fold changes matrix
num_genes <- length(fasta_sequences)
num_groups <- 2 
fold_changes <- matrix(1, nrow=num_genes, ncol=num_groups)

# Update fold changes for the genes in the logfc_table
for (i in 1:nrow(logfc_table)) {
  gene_index <- logfc_table$index[i]
  fc <- logfc_table$fc[i]
  print(gene_index)
  print(num_genes)
  if(gene_index <= num_genes) {
    fold_changes[gene_index, 1] <- fc # Case group
    fold_changes[gene_index, 2] <- 1  # Control group (no change)
  }
}


reads_per_transcript <- rep(reads_per_transcript_count, num_genes) 

# Ensure the reads_per_transcript length is correct
if (length(reads_per_transcript) != num_genes) {
  stop("reads_per_transcript length does not match the number of genes")
}

# Simulation
simulate_experiment(
  fasta = fasta_file,
  num_reps = c(num_reps, num_reps), 
  fold_changes = fold_changes,
  reads_per_transcript = reads_per_transcript,
  outdir = output_dir,
  paired = TRUE
)

# Verify the simulation
simulated_files <- list.files(output_dir, full.names = TRUE)
print(simulated_files)
