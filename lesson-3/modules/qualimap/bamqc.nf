process QUALIMAP_BAMQC {
    conda "bioconda::qualimap=2.2.2d"
    container "biocontainers/qualimap:2.2.2d--1"
    publishDir "${params.outdir}/qualimap", mode: 'copy'

    input:
    tuple path(sorted_bam), path(bai)

    output:
    path "qualimap_results", emit: qc_results

    script:
    """
    qualimap bamqc -bam ${sorted_bam} -outdir qualimap_results -nt ${task.cpus}
    """
}