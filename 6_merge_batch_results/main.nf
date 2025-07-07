#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.batch_dir = "${PWD}/batch_*"
params.output_dir = "${PWD}/merged_results"
params.metadata_dir = "${PWD}/metadata"

// Process to merge expression matrices from all batches
process MERGE_EXPRESSION_MATRICES {
    tag 'merge_expression'
    container 'felixlohmeier/pandas:1.3.3'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        path batch_files

    output:
        path "merged_tpm.csv", emit: tpm
        path "merged_log_tpm.csv", emit: log_tpm
        path "merged_counts.csv", emit: counts
        path "merge_summary.txt", emit: summary

    script:
    """
    python3 ${projectDir}/bin/merge_batch_results.py \
        --batch-dirs ${batch_files} \
        --output-prefix merged
    """
}

// Process to merge samplesheets from all batches
process MERGE_SAMPLESHEETS {
    tag 'merge_samplesheets'
    container 'felixlohmeier/pandas:1.3.3'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        path batch_files

    output:
        path "merged_samplesheet.csv", emit: samplesheet

    script:
    """
    python3 - << 'EOF'
import pandas as pd
import os
import glob

# Find all samplesheet files
samplesheets = []
batch_dirs = "${batch_files}".split()

for batch_dir in batch_dirs:
    samplesheet_path = os.path.join(batch_dir, 'samplesheet.filtered.csv')
    if os.path.exists(samplesheet_path):
        df = pd.read_csv(samplesheet_path)
        df['batch'] = os.path.basename(batch_dir)
        samplesheets.append(df)

if samplesheets:
    merged_df = pd.concat(samplesheets, ignore_index=True)
    merged_df.to_csv('merged_samplesheet.csv', index=False)
    print(f"Merged {len(samplesheets)} samplesheets with {len(merged_df)} total samples")
else:
    print("No samplesheets found to merge")
    # Create empty file
    pd.DataFrame().to_csv('merged_samplesheet.csv', index=False)
EOF
    """
}

// Process to generate normalized log TPM
process NORMALIZE_LOG_TPM {
    tag 'normalize_log_tpm'
    container 'felixlohmeier/pandas:1.3.3'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        path log_tpm

    output:
        path "merged_log_tpm_normalized.csv", emit: normalized

    script:
    """
    python3 - << 'EOF'
import pandas as pd
import numpy as np

# Read log TPM matrix
log_tpm_df = pd.read_csv("${log_tpm}")

# Separate gene column and expression data
gene_col = log_tpm_df.columns[0]
genes = log_tpm_df[gene_col]
expression_data = log_tpm_df.drop(columns=[gene_col])

# Calculate column means
col_means = expression_data.mean()
global_mean = col_means.mean()

# Normalize by subtracting column mean and adding global mean
normalized_data = expression_data.subtract(col_means, axis=1).add(global_mean)

# Combine genes and normalized data
result_df = pd.concat([genes, normalized_data], axis=1)
result_df.to_csv('merged_log_tpm_normalized.csv', index=False)

print(f"Normalized {len(result_df)} genes across {len(expression_data.columns)} samples")
print(f"Global mean: {global_mean:.4f}")
EOF
    """
}

// Process to copy metadata
process COPY_METADATA {
    tag 'copy_metadata'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        path metadata_file

    output:
        path "*.{json,csv}", emit: metadata

    script:
    """
    # Copy all metadata files
    cp ${metadata_file} .
    """
}

// Main workflow
workflow {
    // Find all batch directories
    batch_dirs = Channel
        .fromPath("${params.batch_dir}", type: 'dir', checkIfExists: true)
        .collect()

    // Merge expression matrices
    MERGE_EXPRESSION_MATRICES(batch_dirs)

    // Merge samplesheets
    MERGE_SAMPLESHEETS(batch_dirs)

    // Normalize log TPM
    NORMALIZE_LOG_TPM(MERGE_EXPRESSION_MATRICES.out.log_tpm)

    // Copy metadata files if they exist
    metadata_files = Channel
        .fromPath("${params.metadata_dir}/*.{json,csv}", checkIfExists: false)
        .ifEmpty { file("${params.metadata_dir}/dummy.txt") }

    COPY_METADATA(metadata_files)
}