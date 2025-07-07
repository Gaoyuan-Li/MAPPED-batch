#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.metadata = null
params.sample_ids = null
params.output = null

process CREATE_BATCH_SAMPLESHEET {
    container 'felixlohmeier/pandas:1.3.3'
    
    input:
    path metadata_file
    path sample_ids_file
    
    output:
    path "samplesheet_download.csv", emit: samplesheet
    
    script:
    """
    python3 ${projectDir}/bin/create_batch_samplesheet.py \
        --metadata ${metadata_file} \
        --sample-ids ${sample_ids_file} \
        --output samplesheet_download.csv
    """
}

workflow {
    // Create channels from input parameters
    if (!params.metadata || !params.sample_ids) {
        error "Missing required parameters: metadata and sample_ids must be provided"
    }
    
    // Handle glob pattern for metadata
    metadata_ch = Channel.fromPath(params.metadata, checkIfExists: true)
    sample_ids_ch = Channel.fromPath(params.sample_ids, checkIfExists: true)
    
    // Run the process
    CREATE_BATCH_SAMPLESHEET(metadata_ch, sample_ids_ch)
    
    // Copy output to the specified location
    CREATE_BATCH_SAMPLESHEET.out.samplesheet.subscribe { file ->
        file.copyTo(params.output)
    }
}