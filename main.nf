// Long-read DNA-seq pipeline for Oxford Nanopore Technologies (ONT) data

// CPU path:
// GPU path:

include { refFasta; faidxFor } from './modules/utils/references/references.nf' 
include { PREPARE_SAMPLESHEET_LONG } from './modules/utils/samplesheet/main.nf'


workflow {
    // Parameter validation
    if (!(params.alignment.device in ['cpu', 'gpu'])) {
        error "Invalid alignment device: ${params.alignment.device}. Must be 'cpu' or 'gpu'."
    }

    // Samplesheet preparation
    if (params.samplesheet) {
        samplesheet_ch = channel.fromPath(params.samplesheet, checkIfExists: true)
    }
    else if (params.reads_dir && params.reference_genome) {
        samplesheet_ch = PREPARE_SAMPLESHEET_LONG(params.reads_dir, params.reference_genome)
        samplesheet_ch = PREPARE_SAMPLESHEET_LONG.out.csv
    }
    else {
        error "Either a samplesheet or both reads_dir and reference_genome must be provided."
    }

    // reads_ch first
    reads_ch = samplesheet_ch
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.sample,
                reference: row.reference,
            ]
            tuple(meta, file(row.reads, checkIfExists: true))
        } 

    // Redo basecalling if requested

    // Alignment (+ filtering) step
    
    
}

