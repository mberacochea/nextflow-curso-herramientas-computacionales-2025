process QUALIMAP_BAMQC {
    tag "${meta.id}"
    conda "bioconda::qualimap=2.2.2d"
    container "biocontainers/qualimap:2.2.2d--1"
    publishDir "${params.outdir}/${meta.id}/qualimap", mode: 'copy'

    input:
    tuple val(meta), path(sorted_bam), path(bai)

    output:
    path "${meta.id}_qualimap", emit: qc_results

    script:
    """
    qualimap bamqc -bam ${sorted_bam} -outdir ${meta.id}_qualimap -nt ${task.cpus}
    """
}