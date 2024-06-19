#!/usr/bin/env Rscript

suppressWarnings(library(dplyr))
suppressWarnings(library(polyester))
suppressWarnings(library(Biostrings))

#-------------------------------------------------------------
gene_sets_file <- '../output/gene_sets_with_indices.tsv'
fasta_file <- '../data/genome/gencode.v38.transcripts.fa'
output_dir <- '../output/simulated_rnaseq_data_300'
output_fc_table <- '../output/simulated_rnaseq_data_300/fc_table.tsv'


set_seed <- 123
STN <- 0.2
num_reps <- 3
reads_per_transcript_count <- 300
paired_reads = FALSE

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

# Add FC value to geneset
set.seed(set_seed)

fc_table <- gene_sets %>%
  group_by(gene_set) %>%
  mutate(fc = runif(n(), fc_values[[unique(gene_set)]] - STN, fc_values[[unique(gene_set)]] + STN))

write.table(fc_table, file = output_fc_table, sep = '\t', row.names = FALSE, quote = FALSE)

# Load the FASTA sequences
fasta_sequences <- readDNAStringSet(fasta_file)

# Create a fold changes matrix
num_genes <- length(fasta_sequences)
num_groups <- 2 
fold_changes <- matrix(1, nrow=num_genes, ncol=num_groups)

# Update fold changes for the genes in the fc_table
for (i in 1:nrow(fc_table)) {
  gene_index <- fc_table$index[i]
  fc <- fc_table$fc[i]
  fold_changes[gene_index, 1] <- fc # Case group
  fold_changes[gene_index, 2] <- 1  # Control group (no change)
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
  reads_per_transcript = reads_per_transcript, # None if meanmodel
  #meanmodel=TRUE,
  outdir = output_dir,
  paired = paired_reads
)

# Verify the simulation
simulated_files <- list.files(output_dir, full.names = TRUE)
print(simulated_files)
