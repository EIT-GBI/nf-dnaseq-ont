// Long-read DNA-seq pipeline for Oxford Nanopore Technologies (ONT) data

// CPU path:
// GPU path:

include { DORADO_BASECALLER } from './modules/dorado/basecaller/main.nf'
include { DORADO_TRIM } from './modules/dorado/trim/main.nf'
include { DORADO_ALIGNER } from './modules/dorado/aligner/main.nf'
include { MINIMAP2_ALIGN } from './modules/minimap2/align/main.nf'
include { SAMTOOLS_FASTQ } from './modules/samtools/fastq/main.nf'
include { SAMTOOLS_FILTER } from './modules/samtools/filter/main.nf'
include { SAMTOOLS_MERGE } from './modules/local/samtools/merge/main.nf'
include { SAMTOOLS_FLAGSTAT } from './modules/samtools/flagstat/main.nf'
include { NANOPLOT_NANOPLOT } from './modules/nanoplot/nanoplot/main.nf'
include { MODKIT_PILEUP } from './modules/modkit/pileup/main.nf'
include { MODKIT_MOTIF_SEARCH } from './modules/modkit/motif_search/main.nf'
include { SAMTOOLS_CONSENSUS } from './modules/samtools/consensus/main.nf'
include { BEDTOOLS_BIGWIG } from './modules/bedtools/bigwig/main.nf'


// Samplesheet columns: id, name, ont_r, fasta, fai. One row per (sample,
// barcode) -- a multi-barcode sample repeats id/name/fasta/fai on each row.
// Aliases: sample -> name (or id, when there is no id column), reads -> ont_r,
// reference -> fasta. fai defaults to <fasta>.fai.
// ont_r lives only on the sample record; after the fan-out each meta carries
// one read plus sample/part/n_parts.

def isAbsolutePath(path) {
    path.startsWith('/') || path ==~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\/.*/
}

def resolvePath(value, base, what) {
    def path = "${value ?: ''}".trim()
    if (!path) {
        error "Samplesheet is missing a value for '${what}'."
    }
    if (isAbsolutePath(path)) {
        return file(path, checkIfExists: true)
    }
    if (!base) {
        error "Samplesheet gives a relative path for '${what}' (${path}) but neither params.reference_dir nor a local samplesheet directory is available to resolve it against."
    }
    file("${base.toString().replaceAll('/+$', '')}/${path}", checkIfExists: true)
}

def faidxFor(meta) {
    tuple(meta.fasta, meta.fai)
}

// Sample-level meta: drops the per-barcode fields.
def sampleMeta(meta) {
    meta.subMap(['name', 'fasta', 'fai', 'platform', 'preset', 'n_parts']) + [id: meta.sample]
}


workflow {

    main:
    // Parameter validation
    if (!(params.alignment.tool in ['dorado', 'minimap2'])) {
        error "Invalid alignment tool: ${params.alignment.tool}. Must be 'dorado' or 'minimap2'."
    }

    if (!params.samplesheet) {
        error "params.samplesheet is required."
    }
    samplesheet_ch = channel.fromPath(params.samplesheet, checkIfExists: true)

    // Relative reads resolve against the samplesheet, references against
    // reference_dir.
    def sheet_dir = file(params.samplesheet).parent
    def ref_base  = params.reference_dir ?: sheet_dir

    rows_ch = samplesheet_ch
        .splitCsv(header: true)
        // drop entirely blank rows (trailing ",," lines are common in exported CSVs)
        .filter { row -> row.values().any { v -> v?.toString()?.trim() } }
        .map { row ->
            def id    = row.id?.trim()    ?: row.sample?.trim()
            def reads = row.ont_r?.trim() ?: row.reads?.trim()
            def ref   = row.fasta?.trim() ?: row.reference?.trim()
            if (!id || !reads || !ref) {
                error "Samplesheet needs the columns id,name,ont_r,fasta,fai -- or the older sample,reads,reference -- with one row per barcode. Got: ${row.keySet().join(',')}."
            }
            // resolved now, so a missing index fails in seconds rather than
            // after hours of basecalling
            def fasta = resolvePath(ref, ref_base, "fasta (sample ${id})")
            def fai = row.fai?.trim()
                ? resolvePath(row.fai, ref_base, "fai (sample ${id})")
                : file("${fasta.toUriString()}.fai", checkIfExists: true)

            def info = [
                name:     row.name?.trim() ?: row.sample?.trim() ?: id,
                fasta:    fasta,
                fai:      fai,
                platform: row.platform?.trim() ?: 'ONT',
                preset:   'map-ont',
            ]
            tuple(id, info, resolvePath(reads, sheet_dir, "ont_r (sample ${id})"))
        }

    // One record per sample. Reads are sorted by full path, not basename: in the
    // usual ONT layout every part is called reads.fastq.gz and only the
    // barcodeNN directory tells them apart.
    sample_ch = rows_ch
        .groupTuple()
        .map { id, infos, reads ->
            if (infos.any { info -> info != infos[0] }) {
                error "Samplesheet rows for sample '${id}' disagree on name/fasta/fai/platform -- every barcode row of a sample must repeat the same values."
            }
            def sorted = reads.sort { a, b -> a.toUriString() <=> b.toUriString() }
            // bare map, not tuple(): a one-element tuple arrives as a List,
            // where meta.ont_r is a Groovy spread and nests the read list
            infos[0] + [id: id, ont_r: sorted, n_parts: sorted.size()]
        }

    // Fan out. Whether a sample needs merging is already known (n_parts), so a
    // multi-part sample's files are <id>_1, <id>_2 ... all the way down and the
    // merged BAM takes the plain <id>. A single-part sample never gets a suffix.
    reads_ch = sample_ch
        .flatMap { meta ->
            meta.ont_r.withIndex().collect { read, i ->
                def part_meta = meta.subMap(['name', 'fasta', 'fai', 'platform', 'preset', 'n_parts']) +
                    [id: meta.n_parts > 1 ? "${meta.id}_${i + 1}" : meta.id, sample: meta.id]
                tuple(part_meta, read)
            }
        }

    // Declared up front so publish: can name them even when a branch is off.
    ch_basecalled = channel.empty()
    ch_trimmed    = channel.empty()
    ch_filtered   = channel.empty()
    ch_consensus  = channel.empty()
    ch_bigwig     = channel.empty()
    ch_bedmethyl  = channel.empty()
    ch_motifs     = channel.empty()

    if (params.basecalling.redo) {
        DORADO_BASECALLER(reads_ch, params.basecalling.model)
        ubam_ch       = DORADO_BASECALLER.out.ubam
        ch_basecalled = DORADO_BASECALLER.out.ubam
    }
    else {
        ubam_ch = reads_ch
    }

    if (params.trim.enabled) {
        DORADO_TRIM(ubam_ch)
        trimmed_ch = DORADO_TRIM.out.ubam
        ch_trimmed = DORADO_TRIM.out.ubam
    }
    else {
        trimmed_ch = ubam_ch
    }

    // Alignment. Both tools emit tuple(meta, bam, bai).
    if (params.alignment.tool == 'dorado') {
        aln_in = trimmed_ch.multiMap { meta, reads ->
            reads: tuple(meta, reads)
            fasta: faidxFor(meta)
        }
        DORADO_ALIGNER(aln_in.reads, aln_in.fasta)
        aligned_ch = DORADO_ALIGNER.out.bam
    }
    else {
        // minimap2 cannot read BAM. ext.args carries -T MM,ML here and -y on
        // the aligner, so modified-base tags survive the round trip.
        SAMTOOLS_FASTQ(trimmed_ch)
        aln_in = SAMTOOLS_FASTQ.out.reads.multiMap { meta, reads ->
            reads: tuple(meta, reads)
            fasta: faidxFor(meta)
        }
        MINIMAP2_ALIGN(aln_in.reads, aln_in.fasta)
        aligned_ch = MINIMAP2_ALIGN.out.bam
    }

    if (params.filter.enabled) {
        SAMTOOLS_FILTER(aligned_ch.map { meta, bam, _bai -> tuple(meta, bam) })
        parts_ch    = SAMTOOLS_FILTER.out.filtered
        ch_filtered = SAMTOOLS_FILTER.out.filtered
    }
    else {
        parts_ch = aligned_ch
    }

    // Merge the parts back: one BAM per sample from here on. Single-part
    // samples skip the merge and only shed their per-part meta fields.
    bam_parts = parts_ch.branch { meta, _bam, _bai ->
        single: meta.n_parts == 1
        multi:  true
    }

    // groupKey releases a sample as soon as its own barcodes are done
    merge_in = bam_parts.multi
        .map { meta, bam, _bai -> tuple(groupKey(meta.sample, meta.n_parts), sampleMeta(meta), bam) }
        .groupTuple()
        .map { _key, metas, bams -> tuple(metas[0], bams.sort { a, b -> a.name <=> b.name }) }

    SAMTOOLS_MERGE(merge_in)

    bam_ch = bam_parts.single
        .map { meta, bam, bai -> tuple(sampleMeta(meta), bam, bai) }
        .mix(SAMTOOLS_MERGE.out.bam)

    // The merged BAM is published alongside the parts it was built from rather
    // than in a folder of its own.
    ch_aligned = aligned_ch
    if (params.filter.enabled) {
        ch_filtered = ch_filtered.mix(SAMTOOLS_MERGE.out.bam)
    }
    else {
        ch_aligned = ch_aligned.mix(SAMTOOLS_MERGE.out.bam)
    }

    // QC 
    SAMTOOLS_FLAGSTAT(bam_ch.map {meta, bam, _bai -> tuple(meta, bam) })

    // Coverage track, post-merge so a multi-barcode sample gets one track.
    if (params.bigwig?.enabled) {
        bigwig_in = bam_ch.multiMap { meta, bam, bai ->
            bam:   tuple(meta, bam, bai)
            fasta: faidxFor(meta)
        }
        BEDTOOLS_BIGWIG(bigwig_in.bam, bigwig_in.fasta)
        ch_bigwig = BEDTOOLS_BIGWIG.out.bigwig
    }

    // NANOPLOT_NANOPLOT(bam_ch)

    // Consensus straight from the pileup, no variant calling. ext.args has -a,
    // so sub-threshold positions become N rather than inheriting the reference.
    if (params.consensus?.enabled) {
        SAMTOOLS_CONSENSUS(bam_ch.map { meta, bam, _bai -> tuple(meta, bam) })
        ch_consensus = SAMTOOLS_CONSENSUS.out.fasta
    }

    // Needs basecalling.modified_bases set, so the BAM carries MM/ML tags.
    if (params.methylation?.enabled) {
        pileup_in = bam_ch.multiMap { meta, bam, bai ->
            bam:   tuple(meta, bam, bai)
            fasta: faidxFor(meta)
        }
        MODKIT_PILEUP(pileup_in.bam, pileup_in.fasta)
        ch_bedmethyl = MODKIT_PILEUP.out.bedmethyl

        if (params.methylation?.motif_search) {
            motif_in = MODKIT_PILEUP.out.bedmethyl.multiMap { meta, bed ->
                bed:   tuple(meta, bed)
                fasta: faidxFor(meta)
            }
            MODKIT_MOTIF_SEARCH(motif_in.bed, motif_in.fasta)
            ch_motifs = MODKIT_MOTIF_SEARCH.out.motifs
        }
    }

    publish:
    basecalled  = ch_basecalled
    trimmed    = ch_trimmed
    aligned    = ch_aligned
    filtered   = ch_filtered
    flagstat   = SAMTOOLS_FLAGSTAT.out.flagstat
    consensus  = ch_consensus
    bigwig     = ch_bigwig
    bedmethyl  = ch_bedmethyl
    motifs     = ch_motifs
}

// Published under params.outdir (-output-dir overrides). meta maps in a tuple
// are ignored; every file in it is published.
output {
    basecalled { path 'basecalled'          }
    trimmed    { path 'trimmed'             }
    aligned    { path 'aligned'             }
    filtered   { path 'filtered'            }
    flagstat   { path 'qc/flagstat'         }
    consensus  { path 'consensus'           }
    bigwig     { path 'bigwig'              }
    bedmethyl  { path 'methylation'         }
    motifs     { path 'methylation/motifs'  }
}
