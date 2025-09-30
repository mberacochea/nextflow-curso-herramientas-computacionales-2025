process QUALIMAP_BAMQC {
    cpus 4
    memory '8.GB'
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