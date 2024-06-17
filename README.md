# RNA-seq Simulation Project

This project aims to simulate RNA-seq data with specific differential expression patterns using the `polyester` R package. The simulation process involves extracting gene sets from a JSON file, mapping these genes to their indices in a reference genome FASTA file, and then using these mappings to simulate RNA-seq data with defined fold changes.


## Project Structure

.
├── data  
│   ├── c2.cp.kegg_medicus.v2023.2.Hs.json  
│   └── genome  
│       └── gencode.v38.transcripts.fa  
├── doc   
│   └── environment.yml  
├── indexed_gset.ipynb  
├── License   
├── output  
│   └── gene_sets_with_indices.tsv  
├── README.md  
├── scripts  
│   ├── rnaseq_simulation.R   
│   ├── submit_simulation_EXAMPLE.sh  
│   └── submit_simulation.sh  
└── simulate_RNAseq.Rmd 


## Requirements

- Python 3.x
- Conda
- R with the following packages:
  - dplyr
  - polyester
  - Biostrings


## Requirements

- Python 3.x
- Conda

## Steps and Usage

### Set Up Conda Environment

Create the Conda environment using the `environment.yml` file:

```bash
conda env create -f environment.yml
```

In R session:

BiocManager::install("polyester")


## Steps and Usage

### Step 1: Extract Gene Sets with Indices (Python)

Use the Jupyter notebook `indexed_gset.ipynb` to:
1. Define file paths and target gene sets.
2. Extract gene names from the JSON file.
3. Map gene names to their first occurrence index in the FASTA file.
4. Output a TSV file with the gene set, target gene, and index.

### Step 2: Simulate RNA-seq Data (R)

Use `rnaseq_simulation.Rmd` or the R script `rnaseq_simulation.R` in the `scripts` directory to:
1. Read the TSV file with gene sets and indices.
2. Define fold change values for each gene set.
3. Create a table with fold change values for each gene.
4. Load the FASTA sequences.
5. Create a fold changes matrix.
6. Simulate the RNA-seq experiment with specified replicates and fold changes.

### [OPTIONAL]: Run on Server

Use the `submit_simulation.sh` script in the `scripts` directory to submit the job to the server:
1. From `scripts` directory: submit the job using the command: `sbatch submit_simulation.sh`
