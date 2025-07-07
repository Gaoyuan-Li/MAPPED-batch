#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import glob
import argparse
import sys

def merge_expression_matrices(batch_dirs, output_prefix):
    """
    Merge expression matrices (TPM, log TPM, counts) from multiple batches.
    """
    # Initialize empty dataframes for merging
    all_tpm = None
    all_log_tpm = None
    all_counts = None
    
    # Process each batch directory
    for batch_dir in batch_dirs:
        if not os.path.exists(batch_dir):
            print(f"Warning: Batch directory {batch_dir} not found", file=sys.stderr)
            continue
        
        batch_name = os.path.basename(batch_dir)
        
        # Look for expression matrix files
        # First check in expression_matrices subdirectory
        expr_dir = os.path.join(batch_dir, 'expression_matrices')
        if not os.path.exists(expr_dir):
            # If not found, check in batch_X_expression_matrices format
            expr_dir = batch_dir if 'expression_matrices' in batch_dir else batch_dir
        
        tpm_file = os.path.join(expr_dir, 'tpm.csv')
        log_tpm_file = os.path.join(expr_dir, 'log_tpm.csv')
        counts_file = os.path.join(expr_dir, 'counts.csv')
        
        # Read and merge TPM
        if os.path.exists(tpm_file):
            print(f"Reading TPM from {tpm_file}")
            batch_tpm = pd.read_csv(tpm_file, index_col=0)
            if all_tpm is None:
                all_tpm = batch_tpm
            else:
                # Merge on gene IDs, keeping all genes
                all_tpm = pd.merge(all_tpm, batch_tpm, left_index=True, right_index=True, how='outer')
        
        # Read and merge log TPM
        if os.path.exists(log_tpm_file):
            print(f"Reading log TPM from {log_tpm_file}")
            batch_log_tpm = pd.read_csv(log_tpm_file, index_col=0)
            if all_log_tpm is None:
                all_log_tpm = batch_log_tpm
            else:
                all_log_tpm = pd.merge(all_log_tpm, batch_log_tpm, left_index=True, right_index=True, how='outer')
        
        # Read and merge counts
        if os.path.exists(counts_file):
            print(f"Reading counts from {counts_file}")
            batch_counts = pd.read_csv(counts_file, index_col=0)
            if all_counts is None:
                all_counts = batch_counts
            else:
                all_counts = pd.merge(all_counts, batch_counts, left_index=True, right_index=True, how='outer')
    
    # Fill NaN values with 0 (for genes not present in all batches)
    if all_tpm is not None:
        all_tpm = all_tpm.fillna(0)
    if all_log_tpm is not None:
        all_log_tpm = all_log_tpm.fillna(0)
    if all_counts is not None:
        all_counts = all_counts.fillna(0)
    
    # Save merged expression matrices
    summary = []
    
    if all_tpm is not None:
        output_file = f'{output_prefix}_tpm.csv'
        all_tpm.to_csv(output_file)
        summary.append(f"Merged TPM matrix: {all_tpm.shape[0]} genes x {all_tpm.shape[1]} samples")
        print(f"Saved TPM matrix to {output_file}")
    
    if all_log_tpm is not None:
        output_file = f'{output_prefix}_log_tpm.csv'
        all_log_tpm.to_csv(output_file)
        summary.append(f"Merged log TPM matrix: {all_log_tpm.shape[0]} genes x {all_log_tpm.shape[1]} samples")
        print(f"Saved log TPM matrix to {output_file}")
    
    if all_counts is not None:
        output_file = f'{output_prefix}_counts.csv'
        all_counts.to_csv(output_file)
        summary.append(f"Merged counts matrix: {all_counts.shape[0]} genes x {all_counts.shape[1]} samples")
        print(f"Saved counts matrix to {output_file}")
    
    # Write summary
    with open('merge_summary.txt', 'w') as f:
        f.write('\n'.join(summary))
    
    return all_tpm, all_log_tpm, all_counts

def main():
    parser = argparse.ArgumentParser(description='Merge expression matrices from multiple batches')
    parser.add_argument('--batch-dirs', nargs='+', required=True, help='Batch directories to merge')
    parser.add_argument('--output-prefix', default='merged', help='Prefix for output files')
    
    args = parser.parse_args()
    
    # Merge expression matrices
    merge_expression_matrices(args.batch_dirs, args.output_prefix)

if __name__ == '__main__':
    main()