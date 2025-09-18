#!/usr/bin/env bash
set -euo pipefail

function usage() {
  cat <<EOF
Usage: $0 --organism ORGANISM --outdir OUTDIR --library_layout LIB_LAYOUT --workdir WORKDIR --clean-mode CLEAN_MODE --cpu CPU [--ref-accession REF_ACCESSION] [--max_concurrent_downloads N] [--strain STRAIN]

Options:
  --organism        Organism name (e.g., "Acinetobacter baylyi") - required for metadata download
  --outdir          Output directory for pipeline results
  --workdir         Work directory for Nextflow 'work' files
  --library_layout  Library layout: 'single', 'paired', or 'both'
  --clean-mode      Clean up intermediate files and caches after pipeline completion.
  --cpu             Number of CPUs to allocate per process
  --ref-accession   Optional: specific reference genome accession (e.g., "GCA_008931305.1").
                    If not provided, automatically selects the reference strain for the organism.
  --strain          Optional: filter metadata by strain token in 'ScientificName'.
                    Splits ScientificName on spaces and keeps rows where any token equals
                    or contains the provided string (case-insensitive).
                    Alias: '-strain' also accepted.
  --max_concurrent_downloads  Optional: Maximum number of concurrent downloads (default: 20)
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
STRAIN=""

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
    --strain|-strain)
      STRAIN="$2"
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

run_step2() {
  local batch_label="$1"
  if [[ -n "$batch_label" ]]; then
    echo "=== Step 2: Download FASTQ (${batch_label}) ==="
  else
    echo "=== Step 2: Download FASTQ ==="
  fi
  pushd 2_download_fastq > /dev/null 2>&1
  nextflow run main.nf -work-dir "$WORKDIR" --outdir "$OUTDIR" ${MAX_CONCURRENT_DOWNLOADS:+--max_concurrent_downloads $MAX_CONCURRENT_DOWNLOADS} -resume
  popd > /dev/null 2>&1
}

run_step3() {
  local batch_label="$1"
  if [[ -n "$batch_label" ]]; then
    echo "=== Step 3: Download reference genome (${batch_label}) ==="
  else
    echo "=== Step 3: Download reference genome ==="
  fi
  pushd 3_download_reference_genome > /dev/null 2>&1
  if [[ -n "$REF_ACCESSION" ]]; then
    nextflow run main.nf -work-dir "$WORKDIR" --ref_accession "$REF_ACCESSION" --outdir "$OUTDIR" ${CPU:+--cpu $CPU} -resume
  else
    nextflow run main.nf -work-dir "$WORKDIR" --organism "$ORGANISM" --outdir "$OUTDIR" ${CPU:+--cpu $CPU} -resume
  fi
  popd > /dev/null 2>&1
}

run_step4() {
  local batch_label="$1"
  if [[ -n "$batch_label" ]]; then
    echo "=== Step 4: Generate count/tpm matrix (${batch_label}) ==="
  else
    echo "=== Step 4: Generate count/tpm matrix ==="
  fi
  pushd 4_generate_count_matrix > /dev/null 2>&1
  nextflow run main.nf -work-dir "$WORKDIR" --outdir "$OUTDIR" ${CPU:+--cpu $CPU} -resume
  popd > /dev/null 2>&1
}

# Step 1: Download metadata
echo "=== Step 1: Download metadata ==="
pushd 1_download_metadata_efetch > /dev/null 2>&1
nextflow run main.nf -work-dir "$WORKDIR" --organism "$ORGANISM" --outdir "$OUTDIR" --library_layout "$LIB_LAYOUT" ${STRAIN:+--strain "$STRAIN"} -resume
popd > /dev/null 2>&1

SAMPLE_ID_FILE="$OUTDIR/metadata/sample_id.csv"
if [[ ! -f "$SAMPLE_ID_FILE" ]]; then
  echo "Error: sample_id.csv not found at $SAMPLE_ID_FILE"
  exit 1
fi

ORIGINAL_SAMPLE_ID="$OUTDIR/metadata/sample_id_all.csv"
if [[ ! -f "$ORIGINAL_SAMPLE_ID" ]]; then
  cp "$SAMPLE_ID_FILE" "$ORIGINAL_SAMPLE_ID"
fi

BATCH_SIZE=500
BATCH_DIR="$OUTDIR/metadata/batches"
MANIFEST_PATH="$BATCH_DIR/manifest.json"
mkdir -p "$BATCH_DIR"
rm -f "$BATCH_DIR"/sample_id_batch*.csv "$MANIFEST_PATH"

BATCH_COUNT=$(SAMPLE_ID_FILE="$SAMPLE_ID_FILE" BATCH_DIR="$BATCH_DIR" BATCH_SIZE="$BATCH_SIZE" MANIFEST_PATH="$MANIFEST_PATH" python3 <<'PY'
import csv
import json
import os
import sys

sample_file = os.environ["SAMPLE_ID_FILE"]
batch_dir = os.environ["BATCH_DIR"]
manifest_path = os.environ["MANIFEST_PATH"]
batch_size = int(os.environ.get("BATCH_SIZE", "500"))

with open(sample_file, newline='') as handle:
    reader = csv.reader(handle)
    try:
        header = next(reader)
    except StopIteration:
        header = []
    rows = [row for row in reader if any(cell.strip() for cell in row)]

if not rows:
    with open(manifest_path, 'w') as out_handle:
        json.dump([], out_handle, indent=2)
    print(0)
    sys.exit(0)

manifest = []
for idx, start in enumerate(range(0, len(rows), batch_size), start=1):
    batch_rows = rows[start:start + batch_size]
    batch_name = f"batch_{idx:03d}"
    batch_file = os.path.join(batch_dir, f"sample_id_batch{idx:03d}.csv")

    with open(batch_file, 'w', newline='') as out_handle:
        writer = csv.writer(out_handle)
        if header:
            writer.writerow(header)
        else:
            writer.writerow(['id'])
        writer.writerows(batch_rows)

    samples = []
    for row in batch_rows:
        for cell in row:
            if cell.strip():
                samples.append(cell.strip())
                break
    manifest.append({
        "batch": idx,
        "name": batch_name,
        "file": os.path.basename(batch_file),
        "samples": samples,
    })

with open(manifest_path, 'w') as out_handle:
    json.dump(manifest, out_handle, indent=2)

print(len(manifest))
PY
)

BATCH_COUNT=$(echo "$BATCH_COUNT" | tr -d '[:space:]')
if [[ -z "$BATCH_COUNT" ]]; then
  echo "Error: Failed to determine batch count"
  exit 1
fi

if ! [[ "$BATCH_COUNT" =~ ^[0-9]+$ ]]; then
  echo "Error: Invalid batch count: $BATCH_COUNT"
  exit 1
fi

if [[ "$BATCH_COUNT" -eq 0 ]]; then
  echo "No samples found after metadata download. Exiting."
  exit 0
fi

TOTAL_SAMPLES=$(awk 'NR>1 && NF>0 {count++} END {print count+0}' "$SAMPLE_ID_FILE")
BATCH_OUTPUTS_DIR="$OUTDIR/batch_outputs"
mkdir -p "$BATCH_OUTPUTS_DIR"

echo "Splitting $TOTAL_SAMPLES samples into $BATCH_COUNT batch(es) (max $BATCH_SIZE per batch)"

shopt -s nullglob
mapfile -t BATCH_FILES < <(printf '%s\n' "$BATCH_DIR"/sample_id_batch*.csv | sort)
shopt -u nullglob

if [[ "${#BATCH_FILES[@]}" -ne "$BATCH_COUNT" ]]; then
  echo "Warning: Expected $BATCH_COUNT batch files but found ${#BATCH_FILES[@]}"
fi

REFERENCE_READY=0
BATCH_INDEX=0

for batch_file in "${BATCH_FILES[@]}"; do
  ((BATCH_INDEX++))
  batch_label=$(printf 'Batch %03d' "$BATCH_INDEX")
  batch_dir_name=$(printf 'batch_%03d' "$BATCH_INDEX")
  batch_sample_count=$(awk 'NR>1 && NF>0 {count++} END {print count+0}' "$batch_file")

  echo "=== Processing ${batch_label} (${batch_sample_count} samples) ==="
  cp "$batch_file" "$SAMPLE_ID_FILE"

  # Ensure previous batch artifacts do not block Nextflow publish steps
  if [[ -d "$OUTDIR/samplesheet" ]]; then
    rm -f "$OUTDIR/samplesheet"/samplesheet.csv \
          "$OUTDIR/samplesheet"/samplesheet_download.csv \
          "$OUTDIR/samplesheet"/tmp_samplesheet.csv 2>/dev/null || true
  fi

  run_step2 "$batch_label"

  if [[ "$REFERENCE_READY" -eq 0 ]]; then
    run_step3 "$batch_label"
    REFERENCE_READY=1
  fi

  run_step4 "$batch_label"

  batch_save_dir="$BATCH_OUTPUTS_DIR/$batch_dir_name"
  mkdir -p "$batch_save_dir"

  if [[ -d "$OUTDIR/expression_matrices" ]]; then
    rm -rf "$batch_save_dir/expression_matrices"
    cp -R "$OUTDIR/expression_matrices" "$batch_save_dir/"
  fi

  if [[ -d "$OUTDIR/samplesheet" ]]; then
    rm -rf "$batch_save_dir/samplesheet"
    cp -R "$OUTDIR/samplesheet" "$batch_save_dir/"
  fi

  cp "$batch_file" "$batch_save_dir/sample_id.csv"
  echo "=== Completed ${batch_label} ==="
done

cp "$ORIGINAL_SAMPLE_ID" "$SAMPLE_ID_FILE"

MERGED_EXPRESSION_DIR="$OUTDIR/expression_matrices"
MERGED_SAMPLESHEET_DIR="$OUTDIR/samplesheet"
mkdir -p "$MERGED_EXPRESSION_DIR" "$MERGED_SAMPLESHEET_DIR"

MANIFEST_PATH="$MANIFEST_PATH" BATCH_RESULTS_DIR="$BATCH_OUTPUTS_DIR" OUT_EXPRESSION_DIR="$MERGED_EXPRESSION_DIR" OUT_SAMPLESHEET_DIR="$MERGED_SAMPLESHEET_DIR" python3 <<'PY'
import csv
import json
import os
import sys
from collections import OrderedDict

manifest_path = os.environ['MANIFEST_PATH']
batch_root = os.environ['BATCH_RESULTS_DIR']
expression_out = os.environ['OUT_EXPRESSION_DIR']
samplesheet_out = os.environ['OUT_SAMPLESHEET_DIR']

if not os.path.exists(manifest_path):
    sys.exit(0)

with open(manifest_path) as handle:
    manifest = json.load(handle)

if not manifest:
    sys.exit(0)

manifest = sorted(manifest, key=lambda item: item.get('batch', 0))

def load_matrix(filename):
    gene_ids = None
    matrix_rows = []
    sample_names = []
    for entry in manifest:
        batch_name = entry.get('name')
        if not batch_name:
            continue
        path = os.path.join(batch_root, batch_name, 'expression_matrices', filename)
        if not os.path.exists(path):
            continue
        with open(path, newline='') as fh:
            rows = list(csv.reader(fh))
        if not rows:
            continue
        header = rows[0]
        data_rows = rows[1:]
        current_gene_ids = [row[0] for row in data_rows]
        if gene_ids is None:
            gene_ids = current_gene_ids
            matrix_rows = [[row[0]] + row[1:] for row in data_rows]
        else:
            if current_gene_ids != gene_ids:
                sys.exit(f"Gene order mismatch detected in {path}")
            for idx, row in enumerate(data_rows):
                matrix_rows[idx].extend(row[1:])
        sample_names.extend(header[1:])
    return gene_ids, matrix_rows, sample_names

counts_gene_ids, counts_rows, counts_samples = load_matrix('counts.csv')

if counts_rows:
    if len(counts_samples) != len(set(counts_samples)):
        dupes = sorted({name for name in counts_samples if counts_samples.count(name) > 1})
        sys.exit(f"Duplicate sample columns detected across batches: {dupes}")
    sorted_samples = sorted(counts_samples)
    index_lookup = {name: idx for idx, name in enumerate(counts_samples)}

    def write_matrix(filename, rows):
        path = os.path.join(expression_out, filename)
        with open(path, 'w', newline='') as out_fh:
            writer = csv.writer(out_fh)
            writer.writerow(['GeneID'] + sorted_samples)
            for row in rows:
                values = row[1:]
                ordered = [values[index_lookup[name]] for name in sorted_samples]
                writer.writerow([row[0]] + ordered)

    write_matrix('counts.csv', counts_rows)

    for extra in ['tpm.csv', 'log_tpm.csv']:
        gene_ids, rows, samples = load_matrix(extra)
        if not rows:
            continue
        if gene_ids != counts_gene_ids:
            sys.exit(f"Gene order mismatch detected while merging {extra}")
        if set(samples) != set(counts_samples):
            missing = set(counts_samples) - set(samples)
            extra_samples = set(samples) - set(counts_samples)
            if missing:
                sys.exit(f"Samples missing in {extra}: {sorted(missing)}")
            if extra_samples:
                sys.exit(f"Unexpected samples in {extra}: {sorted(extra_samples)}")
        write_matrix(extra, rows)
else:
    sorted_samples = []


def collect_rows(files, key_field):
    columns = []
    records = OrderedDict()
    for path in files:
        if not os.path.exists(path):
            continue
        with open(path, newline='') as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames:
                for column in reader.fieldnames:
                    if column not in columns:
                        columns.append(column)
            for record in reader:
                if not any((record.get(col, '') or '').strip() for col in reader.fieldnames):
                    continue
                key = (record.get(key_field) or '').strip()
                if not key:
                    continue
                records[key] = record
    return columns, records

# samplesheet.csv
samplesheet_files = [
    os.path.join(batch_root, entry['name'], 'samplesheet', 'samplesheet.csv')
    for entry in manifest
]
columns, records = collect_rows(samplesheet_files, 'sample')
if columns:
    ordered_cols = columns[:]
    if 'sample' in ordered_cols:
        ordered_cols = ['sample'] + [col for col in ordered_cols if col != 'sample']
    else:
        ordered_cols.insert(0, 'sample')
    out_path = os.path.join(samplesheet_out, 'samplesheet.csv')
    with open(out_path, 'w', newline='') as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=ordered_cols)
        writer.writeheader()
        seen = set()
        ordered_keys = sorted_samples if sorted_samples else records.keys()
        for key in ordered_keys:
            key = str(key)
            if key in seen:
                continue
            seen.add(key)
            record = records.get(key)
            if not record:
                continue
            writer.writerow({col: record.get(col, '') for col in ordered_cols})
        for key, record in records.items():
            if key in seen:
                continue
            writer.writerow({col: record.get(col, '') for col in ordered_cols})
elif not records:
    with open(os.path.join(samplesheet_out, 'samplesheet.csv'), 'w', newline='') as out_fh:
        out_fh.write('sample\n')

# samplesheet_download.csv
download_files = [
    os.path.join(batch_root, entry['name'], 'samplesheet', 'samplesheet_download.csv')
    for entry in manifest
]
columns, records = collect_rows(download_files, 'id')
if columns:
    ordered_cols = columns[:]
    if 'id' in ordered_cols:
        ordered_cols = ['id'] + [col for col in ordered_cols if col != 'id']
    else:
        ordered_cols.insert(0, 'id')
    download_order = sorted({sample for entry in manifest for sample in entry.get('samples', [])})
    out_path = os.path.join(samplesheet_out, 'samplesheet_download.csv')
    with open(out_path, 'w', newline='') as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=ordered_cols)
        writer.writeheader()
        seen = set()
        for key in download_order:
            key = str(key)
            if key in seen:
                continue
            seen.add(key)
            record = records.get(key)
            if not record:
                continue
            writer.writerow({col: record.get(col, '') for col in ordered_cols})
        for key, record in records.items():
            if key in seen:
                continue
            writer.writerow({col: record.get(col, '') for col in ordered_cols})
else:
    with open(os.path.join(samplesheet_out, 'samplesheet_download.csv'), 'w', newline='') as out_fh:
        out_fh.write('id\n')

# qc_summary.csv
qc_files = [
    os.path.join(batch_root, entry['name'], 'samplesheet', 'qc_summary.csv')
    for entry in manifest
]
columns, records = collect_rows(qc_files, 'sample')
if columns:
    ordered_cols = columns[:]
    if 'sample' in ordered_cols:
        ordered_cols = ['sample'] + [col for col in ordered_cols if col != 'sample']
    else:
        ordered_cols.insert(0, 'sample')
    out_path = os.path.join(samplesheet_out, 'qc_summary.csv')
    with open(out_path, 'w', newline='') as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=ordered_cols)
        writer.writeheader()
        seen = set()
        ordered_keys = sorted_samples if sorted_samples else records.keys()
        for key in ordered_keys:
            key = str(key)
            if key in seen:
                continue
            seen.add(key)
            record = records.get(key)
            if not record:
                continue
            writer.writerow({col: record.get(col, '') for col in ordered_cols})
        for key, record in records.items():
            if key in seen:
                continue
            writer.writerow({col: record.get(col, '') for col in ordered_cols})
else:
    with open(os.path.join(samplesheet_out, 'qc_summary.csv'), 'w', newline='') as out_fh:
        writer = csv.writer(out_fh)
        writer.writerow(['sample', 'per_base_sequence_quality', 'per_sequence_quality_scores', 'per_base_n_content', 'overall_status'])

# passed_samples.txt
passed_samples = set()
for entry in manifest:
    path = os.path.join(batch_root, entry['name'], 'samplesheet', 'passed_samples.txt')
    if not os.path.exists(path):
        continue
    with open(path) as fh:
        for line in fh:
            value = line.strip()
            if value:
                passed_samples.add(value)

if passed_samples:
    with open(os.path.join(samplesheet_out, 'passed_samples.txt'), 'w') as fh:
        for sample in sorted(passed_samples):
            fh.write(f"{sample}\n")
else:
    open(os.path.join(samplesheet_out, 'passed_samples.txt'), 'w').close()

qc_records = OrderedDict()
qc_csv_path = os.path.join(samplesheet_out, 'qc_summary.csv')
if os.path.exists(qc_csv_path):
    with open(qc_csv_path, newline='') as fh:
        reader = csv.DictReader(fh)
        for record in reader:
            sample = (record.get('sample') or '').strip()
            if not sample:
                continue
            qc_records[sample] = record

total_samples = len(qc_records) if qc_records else len(sorted_samples)
failed_samples = sum(1 for rec in qc_records.values() if (rec.get('overall_status') or '').upper() not in ('', 'PASS'))
with open(os.path.join(samplesheet_out, 'qc_summary.txt'), 'w') as fh:
    fh.write('QC Summary Report\n')
    fh.write('=================\n')
    fh.write(f'Total individual samples processed: {total_samples}\n')
    fh.write(f'Individual samples failed: {failed_samples}\n')
    fh.write(f'Samples passed: {total_samples - failed_samples}\n')

PY

echo "=== Sample Count Summary ==="
if [[ -f "$OUTDIR/samplesheet/samplesheet_download.csv" ]]; then
  download_count=$(awk 'NR>1 && NF>0 {count++} END {print count+0}' "$OUTDIR/samplesheet/samplesheet_download.csv")
  echo "Downloaded experiments (samplesheet_download.csv): $download_count"
else
  echo "samplesheet_download.csv not found"
fi

if [[ -f "$OUTDIR/samplesheet/samplesheet.csv" ]]; then
  filtered_count=$(awk 'NR>1 && NF>0 {count++} END {print count+0}' "$OUTDIR/samplesheet/samplesheet.csv")
  echo "Experiments passing filtration (samplesheet.csv): $filtered_count"
else
  echo "samplesheet.csv not found"
fi
echo "Batches processed: $BATCH_COUNT"
echo "============================="

echo "All batches completed successfully!"

if [[ "$CLEAN_MODE" == "true" ]]; then
  echo "=== Clean mode enabled: cleaning intermediate files ==="

  # Preserve ref_genome folder by moving it to a temporary location
  if [[ -d "$OUTDIR/seqFiles/ref_genome" ]]; then
    echo "Preserving ref_genome folder..."
    mv "$OUTDIR/seqFiles/ref_genome" "$OUTDIR/ref_genome_temp"
  fi

  # Delete everything in OUTDIR except expression_matrices, samplesheet, batch_outputs, and ref_genome_temp
  find "$OUTDIR" -mindepth 1 -maxdepth 1 ! -name expression_matrices ! -name samplesheet ! -name ref_genome_temp ! -name batch_outputs -exec rm -rf {} +

  # Move ref_genome back to the same level as expression_matrices and samplesheet
  if [[ -d "$OUTDIR/ref_genome_temp" ]]; then
    echo "Moving ref_genome to final location..."
    mv "$OUTDIR/ref_genome_temp" "$OUTDIR/ref_genome"
  fi

  # Delete work, .nextflow, and .nextflow.log in module directories
  for sub in 1_download_metadata_efetch 2_download_fastq 3_download_reference_genome 4_generate_count_matrix; do
    rm -rf "$sub/work" "$sub/.nextflow" "$sub/.nextflow.log"
  done

  # Clean the global Nextflow work directory
  rm -rf "$WORKDIR"
fi
