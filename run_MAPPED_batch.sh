#!/usr/bin/env bash
set -euo pipefail

function usage() {
  cat <<EOF
Usage: $0 --organism ORGANISM --outdir OUTDIR --library_layout LIB_LAYOUT --workdir WORKDIR [--clean-mode] --cpu CPU [--ref-accession REF_ACCESSION] [--max_concurrent_downloads N] [--batch_size N]

Options:
  --organism        Organism name (e.g., "Acinetobacter baylyi") - required for metadata download
  --outdir          Output directory for pipeline results
  --workdir         Work directory for Nextflow 'work' files
  --library_layout  Library layout: 'single', 'paired', or 'both'
  --clean-mode      Clean up intermediate files and caches after each batch completion.
  --cpu             Number of CPUs to allocate per process
  --ref-accession   Optional: specific reference genome accession (e.g., "GCA_008931305.1"). 
                    If not provided, automatically selects the reference strain for the organism.
  --max_concurrent_downloads  Optional: Maximum number of concurrent downloads (default: 20)
  --batch_size      Optional: Number of samples per batch (default: 500)
  -h, --help        Show this help message and exit
EOF
}

# Parse arguments
ORGANISM=""
OUTDIR=""
LIB_LAYOUT=""
CLEAN_MODE="false"
CPU=""
WORKDIR=""
REF_ACCESSION=""
MAX_CONCURRENT_DOWNLOADS=""
BATCH_SIZE="500"

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
    --clean-mode)
      CLEAN_MODE="true"
      shift
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
    --batch_size)
      BATCH_SIZE="$2"
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

# Step 1: Download metadata
echo "=== Step 1: Download metadata ==="
pushd 1_download_metadata_efetch > /dev/null 2>&1
nextflow run main.nf -work-dir "$WORKDIR" --organism "$ORGANISM" --outdir "$OUTDIR" --library_layout "$LIB_LAYOUT" -resume
popd > /dev/null 2>&1

# Step 1.5: Create batches
echo "=== Step 1.5: Create batches ==="
pushd 1.5_create_batches > /dev/null 2>&1
nextflow run main.nf -work-dir "$WORKDIR" --outdir "$OUTDIR" --batch_size "$BATCH_SIZE" -resume
popd > /dev/null 2>&1

# Read batch information
if [[ -f "$OUTDIR/batches/batch_info.txt" ]]; then
  echo "=== Batch Information ==="
  cat "$OUTDIR/batches/batch_info.txt"
  echo "========================="
else
  echo "Error: batch_info.txt not found"
  exit 1
fi

# Count number of batch files
BATCH_COUNT=$(ls -1 "$OUTDIR/batches"/batch_*.csv 2>/dev/null | wc -l)
if [[ $BATCH_COUNT -eq 0 ]]; then
  echo "Error: No batch files found"
  exit 1
fi

echo "Found $BATCH_COUNT batches to process"

# Process each batch
BATCH_DIRS=()
for i in $(seq 1 $BATCH_COUNT); do
  BATCH_FILE="$OUTDIR/batches/batch_${i}.csv"
  BATCH_OUTDIR="$OUTDIR/batch_${i}"
  BATCH_WORKDIR="$WORKDIR/batch_${i}"
  
  echo "=== Processing Batch $i of $BATCH_COUNT ==="
  echo "Batch file: $BATCH_FILE"
  echo "Output directory: $BATCH_OUTDIR"
  
  mkdir -p "$BATCH_OUTDIR"
  mkdir -p "$BATCH_WORKDIR"
  
  # Copy metadata to batch directory
  mkdir -p "$BATCH_OUTDIR/metadata"
  cp "$BATCH_FILE" "$BATCH_OUTDIR/metadata/sample_id.csv"
  
  # Also copy the original metadata file for the samples in this batch
  if [[ -f "$OUTDIR/metadata/"*"_metadata.tsv" ]]; then
    cp "$OUTDIR/metadata/"*"_metadata.tsv" "$BATCH_OUTDIR/metadata/"
  fi
  
  # Step 2: Download FASTQ for this batch
  echo "=== Step 2 (Batch $i): Download FASTQ ==="
  pushd 2_download_fastq > /dev/null 2>&1
  nextflow run main.nf -work-dir "$BATCH_WORKDIR" --outdir "$BATCH_OUTDIR" ${MAX_CONCURRENT_DOWNLOADS:+--max_concurrent_downloads $MAX_CONCURRENT_DOWNLOADS} -resume
  popd > /dev/null 2>&1
  
  # Step 3: Download reference genome (only for first batch)
  if [[ $i -eq 1 ]]; then
    echo "=== Step 3 (Batch $i): Download reference genome ==="
    pushd 3_download_reference_genome > /dev/null 2>&1
    if [[ -n "$REF_ACCESSION" ]]; then
      nextflow run main.nf -work-dir "$BATCH_WORKDIR" --ref_accession "$REF_ACCESSION" --outdir "$BATCH_OUTDIR" ${CPU:+--cpu $CPU} -resume
    else
      nextflow run main.nf -work-dir "$BATCH_WORKDIR" --organism "$ORGANISM" --outdir "$BATCH_OUTDIR" ${CPU:+--cpu $CPU} -resume
    fi
    popd > /dev/null 2>&1
  else
    # Copy reference genome from first batch
    echo "=== Copying reference genome from batch 1 to batch $i ==="
    mkdir -p "$BATCH_OUTDIR/seqFiles"
    cp -r "$OUTDIR/batch_1/seqFiles/ref_genome" "$BATCH_OUTDIR/seqFiles/"
  fi
  
  # Step 4: Generate count/tpm matrix for this batch
  echo "=== Step 4 (Batch $i): Generate count/tpm matrix ==="
  pushd 4_generate_count_matrix > /dev/null 2>&1
  nextflow run main.nf -work-dir "$BATCH_WORKDIR" --outdir "$BATCH_OUTDIR" ${CPU:+--cpu $CPU} -resume
  popd > /dev/null 2>&1
  
  # Print sample counts for this batch
  echo "=== Batch $i Sample Count Summary ==="
  if [[ -f "$BATCH_OUTDIR/samplesheet/samplesheet_download.csv" ]]; then
    download_count=$(tail -n +2 "$BATCH_OUTDIR/samplesheet/samplesheet_download.csv" | grep -c '^')
    echo "Downloaded samples (samplesheet_download.csv): $download_count"
  else
    echo "samplesheet_download.csv not found"
  fi
  
  if [[ -f "$BATCH_OUTDIR/samplesheet/samplesheet.csv" ]]; then
    filtered_count=$(tail -n +2 "$BATCH_OUTDIR/samplesheet/samplesheet.csv" | grep -c '^')
    echo "Samples passing filtration (samplesheet.csv): $filtered_count"
  else
    echo "samplesheet.csv not found"
  fi
  echo "============================="
  
  # Add batch directory to list
  BATCH_DIRS+=("$BATCH_OUTDIR")
  
  # Clean mode for this batch if enabled
  if [[ "$CLEAN_MODE" == "true" ]]; then
    echo "=== Clean mode enabled: cleaning intermediate files for batch $i ==="
    
    # Preserve ref_genome, expression_matrices, and samplesheet folders
    if [[ -d "$BATCH_OUTDIR/seqFiles/ref_genome" ]]; then
      echo "Preserving ref_genome folder..."
      mv "$BATCH_OUTDIR/seqFiles/ref_genome" "$BATCH_OUTDIR/ref_genome_temp"
    fi
    
    # Delete everything except expression_matrices, samplesheet, and ref_genome_temp
    find "$BATCH_OUTDIR" -mindepth 1 -maxdepth 1 ! -name expression_matrices ! -name samplesheet ! -name ref_genome_temp -exec rm -rf {} +
    
    # Move ref_genome back
    if [[ -d "$BATCH_OUTDIR/ref_genome_temp" ]]; then
      echo "Moving ref_genome to batch output..."
      mkdir -p "$BATCH_OUTDIR/ref_genome"
      mv "$BATCH_OUTDIR/ref_genome_temp"/* "$BATCH_OUTDIR/ref_genome/"
      rmdir "$BATCH_OUTDIR/ref_genome_temp"
    fi
    
    # Clean the batch work directory
    rm -rf "$BATCH_WORKDIR"
  fi
done

# Step 5: Merge all batches
echo "=== Step 5: Merge all batches ==="
pushd 5_merge_batches > /dev/null 2>&1

# Create batch directories parameter
BATCH_DIRS_PARAM=""
for dir in "${BATCH_DIRS[@]}"; do
  BATCH_DIRS_PARAM="$BATCH_DIRS_PARAM$dir,"
done
BATCH_DIRS_PARAM=${BATCH_DIRS_PARAM%,}  # Remove trailing comma

nextflow run main.nf -work-dir "$WORKDIR" --outdir "$OUTDIR" --batch_dirs "[$BATCH_DIRS_PARAM]" -resume
popd > /dev/null 2>&1

# Final clean up
if [[ "$CLEAN_MODE" == "true" ]]; then
  echo "=== Final cleanup: removing batch directories ==="
  for dir in "${BATCH_DIRS[@]}"; do
    rm -rf "$dir"
  done
  
  # Clean up work directories in each module
  for sub in 1_download_metadata_efetch 1.5_create_batches 2_download_fastq 3_download_reference_genome 4_generate_count_matrix 5_merge_batches; do
    rm -rf "$sub/work" "$sub/.nextflow" "$sub/.nextflow.log*"
  done
  
  # Clean the global Nextflow work directory
  rm -rf "$WORKDIR"
fi

echo "=== Final Sample Count Summary ==="
if [[ -f "$OUTDIR/merged/samplesheet_download.csv" ]]; then
  download_count=$(tail -n +2 "$OUTDIR/merged/samplesheet_download.csv" | grep -c '^')
  echo "Total downloaded samples: $download_count"
else
  echo "merged/samplesheet_download.csv not found"
fi

if [[ -f "$OUTDIR/merged/samplesheet.csv" ]]; then
  filtered_count=$(tail -n +2 "$OUTDIR/merged/samplesheet.csv" | grep -c '^')
  echo "Total samples passing filtration: $filtered_count"
else
  echo "merged/samplesheet.csv not found"
fi
echo "============================="

echo "All steps completed successfully!"
echo "Final output is in: $OUTDIR/merged/"