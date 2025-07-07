#!/usr/bin/env bash
set -euo pipefail

function usage() {
  cat <<EOF
Usage: $0 --organism ORGANISM --outdir OUTDIR --library_layout LIB_LAYOUT --workdir WORKDIR --cpu CPU [--ref-accession REF_ACCESSION] [--max_concurrent_downloads N]

Batch processing version of MAPPED pipeline that processes samples in batches of 500 to manage disk space.

Options:
  --organism        Organism name (e.g., "Acinetobacter baylyi") - required for metadata download
  --outdir          Output directory for pipeline results
  --workdir         Work directory for Nextflow 'work' files
  --library_layout  Library layout: 'single', 'paired', or 'both'
  --cpu             Number of CPUs to allocate per process
  --ref-accession   Optional: specific reference genome accession (e.g., "GCA_008931305.1"). 
                    If not provided, automatically selects the reference strain for the organism.
  --max_concurrent_downloads  Optional: Maximum number of concurrent downloads (default: 20)
  -h, --help        Show this help message and exit
EOF
}

# Parse arguments
ORGANISM=""
OUTDIR=""
LIB_LAYOUT=""
CPU=""
WORKDIR=""
REF_ACCESSION=""
MAX_CONCURRENT_DOWNLOADS=""

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --organism)
      ORGANISM="$2"
      shift 2
      ;;
    --outdir)
      OUTDIR="$2"
      shift 2
      ;;
    --library_layout)
      LIB_LAYOUT="$2"
      shift 2
      ;;
    --workdir)
      WORKDIR="$2"
      shift 2
      ;;
    --cpu)
      CPU="$2"
      shift 2
      ;;
    --ref-accession)
      REF_ACCESSION="$2"
      shift 2
      ;;
    --max_concurrent_downloads)
      MAX_CONCURRENT_DOWNLOADS="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

# Check required arguments
if [[ -z "$ORGANISM" || -z "$OUTDIR" || -z "$LIB_LAYOUT" || -z "$WORKDIR" ]]; then
  echo "Error: Missing required arguments."
  usage
  exit 1
fi

# Validate library_layout parameter
if [[ "$LIB_LAYOUT" != "single" && "$LIB_LAYOUT" != "paired" && "$LIB_LAYOUT" != "both" ]]; then
  echo "Error: Invalid library_layout value: $LIB_LAYOUT"
  echo "Valid values are: single, paired, both"
  exit 1
fi

# Convert OUTDIR to an absolute path and ensure it exists
if [[ "$OUTDIR" != /* ]]; then
  OUTDIR="$(pwd)/$OUTDIR"
fi
mkdir -p "$OUTDIR"

# Convert WORKDIR to an absolute path and ensure it exists
if [[ "$WORKDIR" != /* ]]; then
  WORKDIR="$(pwd)/$WORKDIR"
fi
mkdir -p "$WORKDIR"

# Create batch processing directories
mkdir -p "$OUTDIR/batches"
mkdir -p "$OUTDIR/batch_logs"

echo "=== MAPPED Batch Processing Pipeline ==="
echo "Organism: $ORGANISM"
echo "Output directory: $OUTDIR"
echo "Work directory: $WORKDIR"
echo "Library layout: $LIB_LAYOUT"
echo "CPUs: ${CPU:-default}"
echo "======================================="

# Step 1: Download metadata
echo "=== Step 1: Download metadata ==="
pushd 1_download_metadata_efetch > /dev/null 2>&1
nextflow run main.nf -work-dir "$WORKDIR" --organism "$ORGANISM" --outdir "$OUTDIR" --library_layout "$LIB_LAYOUT" -resume
popd > /dev/null 2>&1

# Check if metadata was downloaded successfully
if [[ ! -f "$OUTDIR/metadata/sample_id.csv" ]]; then
  echo "Error: Metadata download failed. No sample_id.csv found."
  exit 1
fi

# Count total samples
TOTAL_SAMPLES=$(tail -n +2 "$OUTDIR/metadata/sample_id.csv" | wc -l)
echo "Total samples found: $TOTAL_SAMPLES"

if [[ $TOTAL_SAMPLES -eq 0 ]]; then
  echo "Error: No samples found in metadata."
  exit 1
fi

# Step 2: Download reference genome (once for all batches)
echo "=== Step 2: Download reference genome ==="
pushd 2_download_reference_genome > /dev/null 2>&1
if [[ -n "$REF_ACCESSION" ]]; then
  nextflow run main.nf -work-dir "$WORKDIR" --ref_accession "$REF_ACCESSION" --outdir "$OUTDIR" ${CPU:+--cpu $CPU} -resume
else
  nextflow run main.nf -work-dir "$WORKDIR" --organism "$ORGANISM" --outdir "$OUTDIR" ${CPU:+--cpu $CPU} -resume
fi
popd > /dev/null 2>&1

# Check if reference genome was downloaded successfully
if [[ ! -d "$OUTDIR/seqFiles/ref_genome" ]]; then
  echo "Error: Reference genome download failed."
  exit 1
fi

# Calculate number of batches
BATCH_SIZE=500
NUM_BATCHES=$((($TOTAL_SAMPLES + $BATCH_SIZE - 1) / $BATCH_SIZE))
LAST_BATCH_SIZE=$(($TOTAL_SAMPLES % $BATCH_SIZE))

# Adjust batching strategy based on remainder
if [[ $LAST_BATCH_SIZE -ne 0 && $LAST_BATCH_SIZE -lt 250 && $NUM_BATCHES -gt 1 ]]; then
  # Merge last batch with previous if less than 250 samples
  NUM_BATCHES=$(($NUM_BATCHES - 1))
  echo "Adjusting batches: merging last $LAST_BATCH_SIZE samples with previous batch"
fi

echo "Processing $TOTAL_SAMPLES samples in $NUM_BATCHES batches"

# Create batch sample lists
echo "=== Creating batch sample lists ==="
python3 - << EOF
import csv
import os

# Read all sample IDs
with open('$OUTDIR/metadata/sample_id.csv', 'r') as f:
    reader = csv.reader(f)
    header = next(reader)
    sample_ids = [row[0] for row in reader if row]

total_samples = len(sample_ids)
batch_size = 500
num_batches = $NUM_BATCHES

# Create batch directories and sample lists
batch_start = 0
for batch_num in range(1, num_batches + 1):
    batch_dir = f'$OUTDIR/batches/batch_{batch_num}'
    os.makedirs(batch_dir, exist_ok=True)
    
    # Calculate batch end
    if batch_num == num_batches:
        # Last batch gets all remaining samples
        batch_end = total_samples
    else:
        batch_end = min(batch_start + batch_size, total_samples)
    
    # Write batch sample list
    batch_samples = sample_ids[batch_start:batch_end]
    with open(f'{batch_dir}/sample_ids.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['id'])
        for sample_id in batch_samples:
            writer.writerow([sample_id])
    
    print(f"Batch {batch_num}: {len(batch_samples)} samples (indices {batch_start} to {batch_end-1})")
    
    batch_start = batch_end
EOF

# Process each batch
for BATCH_NUM in $(seq 1 $NUM_BATCHES); do
  echo ""
  echo "======================================="
  echo "=== Processing Batch $BATCH_NUM of $NUM_BATCHES ==="
  echo "======================================="
  
  BATCH_DIR="$OUTDIR/batches/batch_${BATCH_NUM}"
  BATCH_WORKDIR="$WORKDIR/batch_${BATCH_NUM}"
  mkdir -p "$BATCH_WORKDIR"
  
  # Step 3: Create batch-specific samplesheet from metadata
  echo "=== Step 3 (Batch $BATCH_NUM): Create batch samplesheet ==="
  
  # Find the metadata file
  METADATA_FILE=$(ls "$OUTDIR/metadata/"*"_metadata.tsv" 2>/dev/null | head -n1)
  if [[ -z "$METADATA_FILE" ]]; then
    echo "Error: No metadata file found in $OUTDIR/metadata/"
    exit 1
  fi
  
  pushd 3_initiate_batch_processing > /dev/null 2>&1
  echo "Debug: Running nextflow with:"
  echo "  metadata: $METADATA_FILE"
  echo "  sample-ids: $BATCH_DIR/sample_ids.csv"
  echo "  output: $BATCH_DIR/samplesheet_download.csv"
  
  nextflow run main.nf \
    -work-dir "$BATCH_WORKDIR" \
    --metadata "$METADATA_FILE" \
    --sample_ids "$BATCH_DIR/sample_ids.csv" \
    --output "$BATCH_DIR/samplesheet_download.csv" \
    -resume
  popd > /dev/null 2>&1
  
  # Step 4: Download FASTQ for this batch
  echo "=== Step 4 (Batch $BATCH_NUM): Download FASTQ ==="
  pushd 4_download_fastq > /dev/null 2>&1
  nextflow run main.nf \
    -work-dir "$BATCH_WORKDIR" \
    --input "$BATCH_DIR/samplesheet_download.csv" \
    --outdir "$BATCH_DIR" \
    ${MAX_CONCURRENT_DOWNLOADS:+--max_concurrent_downloads $MAX_CONCURRENT_DOWNLOADS} \
    -resume
  popd > /dev/null 2>&1
  
  # Copy samplesheet to expected location
  mkdir -p "$BATCH_DIR/samplesheet"
  cp "$BATCH_DIR/samplesheet_download.csv" "$BATCH_DIR/samplesheet/"
  
  # Step 5: Generate count/tpm matrix for this batch
  echo "=== Step 5 (Batch $BATCH_NUM): Generate count/tpm matrix ==="
  
  # Copy reference genome to batch directory to satisfy the pipeline requirements
  mkdir -p "$BATCH_DIR/seqFiles"
  cp -r "$OUTDIR/seqFiles/ref_genome" "$BATCH_DIR/seqFiles/"
  
  pushd 5_generate_count_matrix > /dev/null 2>&1
  nextflow run main.nf \
    -work-dir "$BATCH_WORKDIR" \
    --outdir "$BATCH_DIR" \
    ${CPU:+--cpu $CPU} \
    -resume
  popd > /dev/null 2>&1
  
  # Save batch outputs
  echo "Saving batch $BATCH_NUM outputs..."
  if [[ -d "$BATCH_DIR/expression_matrices" ]]; then
    cp -r "$BATCH_DIR/expression_matrices" "$OUTDIR/batches/batch_${BATCH_NUM}_expression_matrices"
  fi
  if [[ -d "$BATCH_DIR/samplesheet" ]]; then
    cp -r "$BATCH_DIR/samplesheet" "$OUTDIR/batches/batch_${BATCH_NUM}_samplesheet"
  fi
  
  # Clean up batch files to save space
  echo "Cleaning up batch $BATCH_NUM files..."
  # Remove raw and trimmed sequence files
  rm -rf "$BATCH_DIR/seqFiles"
  rm -rf "$BATCH_DIR/fastq"
  rm -rf "$BATCH_DIR/trimmed"
  rm -rf "$BATCH_DIR/salmon"
  rm -rf "$BATCH_DIR/fastqc"
  
  # Clean work directory for this batch
  rm -rf "$BATCH_WORKDIR"
  
  # Save batch log
  if [[ -f 4_generate_count_matrix/.nextflow.log ]]; then
    cp 4_generate_count_matrix/.nextflow.log "$OUTDIR/batch_logs/batch_${BATCH_NUM}.log"
  fi
  
  echo "Batch $BATCH_NUM completed and cleaned."
done

# Step 6: Merge all batch results
echo ""
echo "======================================="
echo "=== Step 6: Merging batch results ==="
echo "======================================="

# Run the merge batch results Nextflow pipeline
cd 6_merge_batch_results
nextflow run main.nf \
  --batch_dir "$OUTDIR/batches/batch_*_expression_matrices" \
  --output_dir "$OUTDIR/merged_results" \
  --metadata_dir "$OUTDIR/metadata" \
  -profile docker

# Check if the pipeline completed successfully
if [ $? -eq 0 ]; then
  echo "Merge batch results pipeline completed successfully"
  
  # Copy merged results to the expected locations
  mkdir -p "$OUTDIR/expression_matrices"
  if [ -f "$OUTDIR/merged_results/merged_tpm.csv" ]; then
    cp "$OUTDIR/merged_results/merged_tpm.csv" "$OUTDIR/expression_matrices/tpm.csv"
  fi
  if [ -f "$OUTDIR/merged_results/merged_log_tpm.csv" ]; then
    cp "$OUTDIR/merged_results/merged_log_tpm.csv" "$OUTDIR/expression_matrices/log_tpm.csv"
  fi
  if [ -f "$OUTDIR/merged_results/merged_counts.csv" ]; then
    cp "$OUTDIR/merged_results/merged_counts.csv" "$OUTDIR/expression_matrices/counts.csv"
  fi
  if [ -f "$OUTDIR/merged_results/merged_log_tpm_normalized.csv" ]; then
    cp "$OUTDIR/merged_results/merged_log_tpm_normalized.csv" "$OUTDIR/expression_matrices/log_tpm_norm.csv"
  fi
  
  # Copy merged samplesheet
  mkdir -p "$OUTDIR/samplesheet"
  if [ -f "$OUTDIR/merged_results/merged_samplesheet.csv" ]; then
    cp "$OUTDIR/merged_results/merged_samplesheet.csv" "$OUTDIR/samplesheet/samplesheet.csv"
  fi
else
  echo "Error: Merge batch results pipeline failed"
  exit 1
fi

cd ..

# Copy reference genome to final location
echo "Copying reference genome to final location..."
cp -r "$OUTDIR/seqFiles/ref_genome" "$OUTDIR/ref_genome"

# Clean up batch directories
echo "Cleaning up batch directories..."
rm -rf "$OUTDIR/batches"
rm -rf "$OUTDIR/seqFiles"

# Clean up Nextflow logs
for sub in 1_download_metadata_efetch 2_download_reference_genome 3_initiate_batch_processing 4_download_fastq 5_generate_count_matrix 6_merge_batch_results; do
  rm -rf "$sub/.nextflow" "$sub/.nextflow.log"*
done

echo ""
echo "======================================="
echo "=== MAPPED batch processing completed! ==="
echo "======================================="
echo "Results saved in: $OUTDIR"
echo "- Expression matrices: $OUTDIR/expression_matrices/"
echo "- Sample sheet: $OUTDIR/samplesheet/"
echo "- Reference genome: $OUTDIR/ref_genome/"
echo "- Metadata: $OUTDIR/metadata/"
echo "- Batch logs: $OUTDIR/batch_logs/"