# MAPPED-batch - Modular Automated Pipeline for Public Expression Data (Batch Version)

MAPPED-batch is a comprehensive Nextflow-based workflow designed to analyze large-scale public RNA-seq data from NCBI SRA in batches. It automates the entire process from metadata retrieval to expression matrix generation, with built-in batch processing capabilities for handling thousands of samples efficiently.

## Overview

MAPPED-batch consists of five integrated modules that work together to process public expression data in batches:

1. **Metadata Download**: Retrieves and formats metadata from NCBI SRA based on organism name
2. **Batch Creation**: Divides samples into manageable batches (default: 500 samples per batch)
3. **FASTQ Download**: Efficiently downloads sequencing data for each batch
4. **Reference Genome Download**: Obtains reference genome sequences and annotations
5. **Expression Quantification**: Performs quality control, trimming, and gene expression quantification per batch
6. **Batch Merging**: Combines results from all batches and performs final normalization

The pipeline is designed to handle large-scale datasets with built-in batch processing, error handling, resume capabilities, and resource optimization.

## Prerequisites

- **[Nextflow](https://www.nextflow.io/)** (version 21.04.0 or later)
- **[Docker](https://www.docker.com/)** (version 20.10 or later)

## Installation

1. Clone the MAPPED repository:
```bash
git clone https://github.com/your-org/MAPPED-batch.git
cd MAPPED-batch
```

2. Ensure the wrapper script is executable:
```bash
chmod +x run_MAPPED_batch.sh
```

3. Verify Nextflow and Docker are installed:
```bash
nextflow -version
docker --version
```

## Quick Start

Process RNA-seq data for an organism using the default reference genome:

```bash
./run_MAPPED_batch.sh \
    --organism "Escherichia coli" \
    --outdir ./results \
    --workdir ./work \
    --library_layout paired \
    --cpu 48
```

## Usage

### Basic Usage

The `run_MAPPED_batch.sh` wrapper script orchestrates all pipeline modules with batch processing:

```bash
./run_MAPPED.sh [OPTIONS]
```

### Using a Specific Reference Genome

To use a specific genome assembly instead of the default reference strain:

```bash
./run_MAPPED_batch.sh \
    --organism "Streptomyces coelicolor" \
    --ref-accession GCA_008931305.1 \
    --outdir ./results \
    --workdir ./work \
    --library_layout paired \
    --cpu 24
```

### Clean Mode

To automatically clean up intermediate files after successful completion:

```bash
./run_MAPPED_batch.sh \
    --organism "Pseudomonas putida" \
    --outdir ./results \
    --workdir ./work \
    --library_layout paired \
    --cpu 16 \
    --clean-mode
```

## Pipeline Modules

### 1. Download Metadata (Module 1)
- Queries NCBI SRA for RNA-seq experiments matching the specified organism
- Filters samples based on library layout (single-end, paired-end, or both)
- Generates formatted metadata files for downstream processing

### 2. Create Batches (Module 1.5)
- Divides all samples into batches of configurable size (default: 500)
- Intelligently handles the last batch:
  - If remaining samples < 250: merges with previous batch
  - If remaining samples ≥ 250: creates separate batch

### 3. Download FASTQ (Module 2)
- Downloads raw sequencing data for each batch independently
- Validates downloaded files
- Creates a samplesheet for downstream analysis

### 4. Download Reference Genome (Module 3)
- Downloads reference genome assemblies from NCBI (only for first batch)
- Retrieves genome sequence (FASTA), annotations (GFF), and protein sequences (FAA)
- Supports two modes:
  - **Default mode**: Automatically selects the largest reference genome for the organism
  - **Accession mode**: Downloads a specific genome assembly using its accession number

### 5. Generate Count Matrix (Module 4)
- Performs quality control on raw reads (FastQC) per batch
- Trims adapters and low-quality bases (TrimGalore)
- Quantifies gene expression using Salmon
- Generates expression matrices (TPM, log TPM, and raw counts)
- Creates comprehensive quality reports (MultiQC)
- **Note**: Normalization is skipped at batch level

### 6. Merge Batches (Module 5)
- Combines samplesheets from all batches
- Merges expression matrices (TPM, log TPM, counts)
- Performs final normalization on the merged log TPM data
- Copies reference genome to final output

## Parameters

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--organism` | Full taxonomic name of the target organism | `"Escherichia coli"` |
| `--outdir` | Output directory for all results | `/path/to/results` |
| `--workdir` | Nextflow work directory for temporary files | `/path/to/work` |
| `--library_layout` | Sequencing library type: `single`, `paired`, or `both` | `paired` |

### Optional Parameters

| Parameter | Description | Default | Example |
|-----------|-------------|---------|---------|
| `--cpu` | Number of CPUs to allocate per process | System dependent | `16` |
| `--ref-accession` | Specific reference genome accession | Auto-selected | `GCA_008931305.1` |
| `--max_concurrent_downloads` | Maximum number of concurrent FASTQ downloads | `20` | `10` |
| `--batch_size` | Number of samples per batch | `500` | `1000` |
| `--clean-mode` | Remove intermediate files after each batch | `false` | (flag) |
| `-h, --help` | Display help message | - | (flag) |

## Output Structure

The pipeline creates a well-organized output directory structure with batch processing artifacts:

```
${outdir}/
├── metadata/                    # Downloaded and formatted metadata
│   ├── *_metadata.tsv           # Complete metadata for all samples
│   └── sample_id.csv            # List of SRA accessions
├── batches/                     # Batch information
│   ├── batch_1.csv              # Sample IDs for batch 1
│   ├── batch_2.csv              # Sample IDs for batch 2
│   ├── ...                      # Additional batch files
│   └── batch_info.txt           # Summary of batch sizes
├── batch_1/                     # Results for batch 1 (if not cleaned)
│   ├── samplesheet/
│   ├── expression_matrices/
│   └── ...
├── batch_2/                     # Results for batch 2 (if not cleaned)
│   └── ...
├── ...                          # Additional batch directories
└── merged/                      # Final merged results
    ├── samplesheet/             # Combined sample information
    │   ├── samplesheet_download.csv # All available samples from NCBI
    │   └── samplesheet.csv      # Samples that passed QC
    ├── ref_genome/              # Reference genome files
    │   ├── *.fna                # Genome sequence (FASTA)
    │   ├── *.gff                # Gene annotations (GFF3)
    │   ├── *.faa                # Protein sequences (FASTA)
    │   └── datasets_summary.json
    └── expression_matrices/     # Final merged expression data
        ├── counts.csv           # Raw count matrix
        ├── tpm.csv              # TPM normalized matrix
        ├── log_tpm.csv          # Log-transformed TPM
        └── log_tpm_norm.csv     # Log-transformed normalized TPM
```

### Clean Mode Output

When using `--clean-mode`, intermediate batch files are cleaned after processing:
```
${outdir}/
├── metadata/               # Original metadata
├── batches/                # Batch information files
└── merged/                 # Final merged results
    ├── expression_matrices/     # Final expression matrices
    ├── samplesheet/            # Combined sample metadata
    └── ref_genome/             # Reference genome files
```

## Batch Processing Details

### Batch Size Logic
- Default batch size: 500 samples
- Last batch handling:
  - If remaining samples < 250: merge with previous batch
  - If remaining samples ≥ 250: create separate batch
  - Example: 1200 samples → batch 1 (500), batch 2 (700)
  - Example: 1400 samples → batch 1 (500), batch 2 (500), batch 3 (400)

### Clean Mode Behavior
When `--clean-mode` is enabled:
- After each batch completes: removes all intermediate files except expression matrices, samplesheets, and reference genome
- After final merge: removes individual batch directories
- Significantly reduces disk usage for large datasets