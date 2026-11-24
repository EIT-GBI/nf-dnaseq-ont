// Merge the per-barcode BAMs of one sample. Inputs must be coordinate-sorted
// (both aligners sort). samtools merge keeps aux tags, so MM/ML survive.
// Vendored because nf-mod-samtools has no merge module yet.
// TODO: upstream into eit-gbi/nf-mod-samtools.
//
// KNOWN GAP: the aligners build @RG from meta.id, which here is the per-part
// id, so the merged BAM has one @RG per part with SM:<sample>_N.
// Nothing downstream in this pipeline groups by SM, but a tool that does would
// see two samples. Fix belongs in nf-mod-minimap2/nf-mod-dorado: SM from
// `meta.sample ?: meta.id`, barcode in LB.

process SAMTOOLS_MERGE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam

    script:
    def args = task.ext.args ?: ''
    """
    samtools merge \\
        ${args} \\
        -@ ${task.cpus} \\
        --write-index \\
        -o ${meta.id}.bam##idx##${meta.id}.bam.bai \\
        ${bams}
    """

    stub:
    """
    touch ${meta.id}.bam ${meta.id}.bam.bai
    """
}
