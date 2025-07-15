#!/usr/bin/env nextflow

params.outdir = null
params.batch_dirs = null  // Will be passed as a list of batch directories

process MERGE_BATCH_DATA {
    tag 'merge_batches'
    container 'felixlohmeier/pandas:1.3.3'
    publishDir "${params.outdir}/merged", mode: 'copy', overwrite: true
    
    input:
        path batch_dirs
    
    output:
        path 'samplesheet.csv', emit: samplesheet
        path 'samplesheet_download.csv', emit: samplesheet_download
        path 'tpm.csv', emit: tpm
        path 'log_tpm.csv', emit: log_tpm
        path 'counts.csv', emit: counts
    
    script:
    """
    python3 - << 'EOF'
import pandas as pd
import os
import glob

print("=== MERGE_BATCH_DATA ===")

batch_dirs = '${batch_dirs}'.split()
print(f"Merging data from {len(batch_dirs)} batches...")

# Initialize empty dataframes
merged_samplesheet = []
merged_samplesheet_download = []
merged_tpm = None
merged_log_tpm = None
merged_counts = None

for batch_dir in batch_dirs:
    print(f"\\nProcessing batch: {batch_dir}")
    
    # Read samplesheets
    samplesheet_path = os.path.join(batch_dir, 'samplesheet', 'samplesheet.csv')
    if os.path.exists(samplesheet_path):
        df = pd.read_csv(samplesheet_path)
        merged_samplesheet.append(df)
        print(f"  - Found {len(df)} samples in samplesheet.csv")
    
    samplesheet_download_path = os.path.join(batch_dir, 'samplesheet', 'samplesheet_download.csv')
    if os.path.exists(samplesheet_download_path):
        df = pd.read_csv(samplesheet_download_path)
        merged_samplesheet_download.append(df)
        print(f"  - Found {len(df)} samples in samplesheet_download.csv")
    
    # Read expression matrices
    tpm_path = os.path.join(batch_dir, 'expression_matrices', 'tpm.csv')
    if os.path.exists(tpm_path):
        df = pd.read_csv(tpm_path, index_col=0)
        if merged_tpm is None:
            merged_tpm = df
        else:
            # Merge columns (samples)
            merged_tpm = pd.concat([merged_tpm, df], axis=1)
        print(f"  - TPM matrix: {df.shape[1]} samples, {df.shape[0]} genes")
    
    log_tpm_path = os.path.join(batch_dir, 'expression_matrices', 'log_tpm.csv')
    if os.path.exists(log_tpm_path):
        df = pd.read_csv(log_tpm_path, index_col=0)
        if merged_log_tpm is None:
            merged_log_tpm = df
        else:
            merged_log_tpm = pd.concat([merged_log_tpm, df], axis=1)
        print(f"  - Log TPM matrix: {df.shape[1]} samples")
    
    counts_path = os.path.join(batch_dir, 'expression_matrices', 'counts.csv')
    if os.path.exists(counts_path):
        df = pd.read_csv(counts_path, index_col=0)
        if merged_counts is None:
            merged_counts = df
        else:
            merged_counts = pd.concat([merged_counts, df], axis=1)
        print(f"  - Counts matrix: {df.shape[1]} samples")

# Combine samplesheets
print("\\nMerging samplesheets...")
if merged_samplesheet:
    final_samplesheet = pd.concat(merged_samplesheet, ignore_index=True)
    final_samplesheet.to_csv('samplesheet.csv', index=False)
    print(f"Final samplesheet: {len(final_samplesheet)} samples")

if merged_samplesheet_download:
    final_samplesheet_download = pd.concat(merged_samplesheet_download, ignore_index=True)
    final_samplesheet_download.to_csv('samplesheet_download.csv', index=False)
    print(f"Final samplesheet_download: {len(final_samplesheet_download)} samples")

# Save merged expression matrices
print("\\nSaving merged expression matrices...")
if merged_tpm is not None:
    merged_tpm.to_csv('tpm.csv')
    print(f"TPM matrix: {merged_tpm.shape}")

if merged_log_tpm is not None:
    merged_log_tpm.to_csv('log_tpm.csv')
    print(f"Log TPM matrix: {merged_log_tpm.shape}")

if merged_counts is not None:
    merged_counts.to_csv('counts.csv')
    print(f"Counts matrix: {merged_counts.shape}")

print("\\nMerge complete!")
EOF
    """
}

process NORMALIZE_MERGED_LOG_TPM {
    tag 'normalize_merged_log_tpm'
    container 'felixlohmeier/pandas:1.3.3'
    publishDir "${params.outdir}/merged", mode: 'copy', overwrite: true
    
    input:
        path log_tpm_csv
    
    output:
        path 'log_tpm_norm.csv'
    
    script:
    """
    python3 - << 'EOF'
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def normalize_csv(csv_path_in: str, csv_path_out: str):
    # Normalize a CSV file by subtracting the row-wise mean of all
    # columns to set each gene's "average" expression to zero.
    df = pd.read_csv(csv_path_in, index_col=0)
    print(f"Loaded data with shape: {df.shape}")
    
    # Calculate row means
    row_means = df.mean(axis=1)
    
    # Subtract row means from each value in the row
    df_normalized = df.subtract(row_means, axis=0)
    
    # Save to output
    df_normalized.to_csv(csv_path_out)
    print(f"Normalized data saved to {csv_path_out}")
    
    # Print some statistics
    print(f"Mean of row means: {row_means.mean():.4f}")
    print(f"Std of row means: {row_means.std():.4f}")
    print(f"Mean of normalized data: {df_normalized.mean().mean():.4e} (should be ~0)")

normalize_csv('${log_tpm_csv}', 'log_tpm_norm.csv')
EOF
    """
}

process COPY_REF_GENOME {
    tag 'copy_ref_genome'
    publishDir "${params.outdir}/merged", mode: 'copy', overwrite: true
    
    input:
        path ref_genome_dir
    
    output:
        path 'ref_genome/*'
    
    script:
    """
    mkdir -p ref_genome
    cp -r ${ref_genome_dir}/* ref_genome/
    """
}

workflow {
    if (!params.outdir || !params.batch_dirs) {
        error "Please provide --outdir and --batch_dirs parameters"
    }
    
    // Convert batch directories to channel
    batch_dirs_ch = Channel.fromList(params.batch_dirs)
        .collect()
    
    // Merge batch data
    merged_data = MERGE_BATCH_DATA(batch_dirs_ch)
    
    // Normalize the merged log TPM
    NORMALIZE_MERGED_LOG_TPM(merged_data.log_tpm)
    
    // Copy reference genome from first batch
    ref_genome_ch = Channel.fromPath("${params.batch_dirs[0]}/ref_genome")
    COPY_REF_GENOME(ref_genome_ch)
}

// Clean up log files on completion
workflow.onComplete {
    def logPattern = ~/\.nextflow\.log\.\d+/  
    new File('.').listFiles().findAll { it.name ==~ logPattern }.each { it.delete() }
}