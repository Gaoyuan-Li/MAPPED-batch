#!/usr/bin/env nextflow

def runStartTime = System.currentTimeMillis()

params.organism = null
params.outdir = null

process FETCH_METADATA {

    publishDir 'tmp', mode: 'copy'

    container 'quay.io/biocontainers/entrez-direct:22.4--he881be0_0'

    input:
        val organism
    output:

        path 'tmp_metadata.tsv'

    script:

        def query = '"' + organism + '"[Organism] AND "rna seq"[Strategy] AND "transcriptomic"[Source]'
        """
        esearch -db sra -query '${query}' | efetch -db sra -format runinfo > tmp_metadata.tsv
        """
}

process FORMAT_METADATA {

    publishDir "${params.outdir}/metadata", mode: 'copy'   // ⬅ copy results into metadata subfolder of outdir

    container 'felixlohmeier/pandas:1.3.3'

    input:
        path raw_tsv
        path clean_script
        val  organism
        val  library_layout

    output:
        path "*_metadata.tsv"                 // ⬅ any metadata file
        path "sample_id.csv"                  // ⬅ sample IDs for downstream use

    script:
        def safe_name = organism.replaceAll(/\s+/, '_')
        def outfile   = "${safe_name}_metadata.tsv"

        """
        python3 ${clean_script} -i ${raw_tsv} -o ${outfile} -l ${library_layout}
        """
}

workflow {
    if ( !params.organism || !params.outdir ) {
        error "You must provide both --organism and --outdir parameters."
    }

    raw_metadata = FETCH_METADATA( params.organism )

    clean_script = file( 'bin/clean_metadata_file.py' )

    ( cleaned_metadata, sample_ids ) = FORMAT_METADATA(
        raw_metadata,
        clean_script,
        params.organism,
        params.library_layout
    )

}

// Add an onComplete event handler to always delete rotated Nextflow log files
workflow.onComplete {
    def logPattern = ~/\.nextflow\.log\.\d+/  
    new File('.').listFiles().findAll { it.name ==~ logPattern }.each { it.delete() }
    // delete the tmp directory
    def tmpDir = new File('tmp')
    if (tmpDir.exists()) {
        tmpDir.deleteDir()
    }
}