// Long-read DNA-seq pipeline for Oxford Nanopore Technologies (ONT) data

// CPU path:
// GPU path:

include { refFasta; faidxFor } from './modules/utils/references/references.nf' 
include { PREPARE_SAMPLESHEET_LONG } from './modules/utils/samplesheet/main.nf'

include { DORADO_BASECALLER } from './modules/dorado/basecaller/main.nf'
include { DORADO_TRIM } from './modules/dorado/trim/main.nf'
include { DORADO_ALIGNER } from './modules/dorado/aligner/main.nf'
include { SAMTOOLS_FILTER } from './modules/samtools/filter/main.nf'
include { SAMTOOLS_FLAGSTAT } from './modules/samtools/flagstat/main.nf'
include { NANOPLOT_NANOPLOT } from './modules/nanoplot/nanoplot/main.nf'
include { MODKIT_PILEUP } from './modules/modkit/pileup/main.nf'
include { MODKIT_MOTIF_SEARCH } from './modules/modkit/motif_search/main.nf'


workflow {
    // Parameter validation
    if (!(params.alignment.device in ['cpu', 'gpu'])) {
        error "Invalid alignment device: ${params.alignment.device}. Must be 'cpu' or 'gpu'."
    }
    if (params.alignment.tool != 'dorado') {
        error "Invalid alignment tool: ${params.alignment.tool}. Only 'dorado' is implemented."
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
        // drop entirely blank rows (trailing ",,"" lines are common in exported CSVs)
        .filter { row -> row.values().any { it?.toString()?.trim() } }
        .map { row ->
            def meta = [
                id: row.sample,
                reference: row.reference,
            ]
            // resolve the reference + .fai now, so a missing index fails in
            // seconds rather than after basecalling has run for hours
            faidxFor(meta)
            tuple(meta, file(row.reads, checkIfExists: true))
        } 

    // Redo basecalling if requested
    if (params.basecalling.redo) {
        DORADO_BASECALLER(reads_ch)
        ubam_ch = DORADO_BASECALLER.out.ubam
    }
    else {
        ubam_ch = reads_ch
    }

    // Alignment (+ filtering) step
    if (params.alignment.tool == 'dorado') {
        if (params.trim.enabled) {
            DORADO_TRIM(ubam_ch)
            trimmed_ch = DORADO_TRIM.out.ubam
        }
        else {
            trimmed_ch = ubam_ch
        }

        aln_in = trimmed_ch.multiMap { meta, reads ->
            reads: tuple(meta, reads)
            fasta: faidxFor(meta)
        }
        DORADO_ALIGNER(aln_in.reads, aln_in.fasta)

        // Do filtering if requested
        if (params.filter.enabled) {
            SAMTOOLS_FILTER(DORADO_ALIGNER.out.bam.map {meta, bam, _bai -> tuple(meta, bam) })
            bam_ch = SAMTOOLS_FILTER.out.filtered
        }
        else {
            bam_ch = DORADO_ALIGNER.out.bam
        }
    }

    // QC 
    SAMTOOLS_FLAGSTAT(bam_ch.map {meta, bam, _bai -> tuple(meta, bam) })

    NANOPLOT_NANOPLOT(bam_ch)

    // Methylation: per-site pileup, then de novo motif discovery.
    // Requires basecalling.modified_bases to have been set, so the BAM
    // carries MM/ML tags.
    if (params.methylation?.enabled) {
        pileup_in = bam_ch.multiMap { meta, bam, bai ->
            bam:   tuple(meta, bam, bai)
            fasta: faidxFor(meta)
        }
        MODKIT_PILEUP(pileup_in.bam, pileup_in.fasta)

        if (params.methylation?.motif_search) {
            motif_in = MODKIT_PILEUP.out.bedmethyl.multiMap { meta, bed ->
                bed:   tuple(meta, bed)
                fasta: faidxFor(meta)
            }
            MODKIT_MOTIF_SEARCH(motif_in.bed, motif_in.fasta)
        }
    }

}

