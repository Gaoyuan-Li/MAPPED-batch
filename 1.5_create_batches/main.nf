#!/usr/bin/env nextflow

params.outdir = null
params.batch_size = 500

process CREATE_BATCHES {
    tag 'create_batches'
    container 'felixlohmeier/pandas:1.3.3'
    publishDir "${params.outdir}/batches", mode: 'copy'

    input:
        path sample_ids

    output:
        path "batch_*.csv", emit: batch_files
        path "batch_info.txt", emit: batch_info

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import math

    # Read sample IDs
    samples_df = pd.read_csv('${sample_ids}')
    total_samples = len(samples_df)
    batch_size = ${params.batch_size}
    
    # Calculate number of batches
    num_batches = math.ceil(total_samples / batch_size)
    
    # Check if we need to merge the last two batches
    last_batch_size = total_samples % batch_size
    if last_batch_size > 0 and last_batch_size < 250 and num_batches > 1:
        # Merge last two batches
        num_batches -= 1
    
    # Create batches
    batch_info = []
    for i in range(num_batches):
        start_idx = i * batch_size
        if i == num_batches - 1:
            # Last batch takes all remaining samples
            end_idx = total_samples
        else:
            end_idx = min((i + 1) * batch_size, total_samples)
        
        batch_df = samples_df.iloc[start_idx:end_idx]
        batch_filename = f'batch_{i+1}.csv'
        batch_df.to_csv(batch_filename, index=False)
        
        batch_info.append(f"Batch {i+1}: {len(batch_df)} samples ({start_idx+1}-{end_idx})")
    
    # Write batch information
    with open('batch_info.txt', 'w') as f:
        f.write(f"Total samples: {total_samples}\\n")
        f.write(f"Number of batches: {num_batches}\\n")
        f.write("\\n".join(batch_info))
    """
}

workflow {
    if (!params.outdir) {
        error "Please provide --outdir parameter"
    }
    
    sample_ids = Channel.fromPath("${params.outdir}/metadata/sample_id.csv")
    CREATE_BATCHES(sample_ids)
}