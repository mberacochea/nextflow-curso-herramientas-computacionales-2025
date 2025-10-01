process SAMTOOLS_SORT {
    tag "${meta.id}"
    cpus 2
    memory '1.GB'
    conda "bioconda::samtools=1.17"
    container "quay.io/biocontainers/samtools:1.17--hd87286a_2"
    publishDir "${params.outdir}/${meta.id}/bam", mode: 'copy', pattern: "*.sorted.bam"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), emit: sorted_bam

    script:
    """
    samtools sort -@ ${task.cpus} ${bam} -o ${meta.id}.sorted.bam
    """
}