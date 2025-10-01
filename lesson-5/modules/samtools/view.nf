process SAMTOOLS_VIEW {
    tag "${meta.id}"
    conda "bioconda::samtools=1.17"
    container "biocontainers/samtools:1.17--hd87286a_2"
    publishDir "${params.outdir}/${meta.id}/bam", mode: 'copy', pattern: "*.bam"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam

    script:
    """
    samtools view -@ ${task.cpus} -bS ${sam} > ${meta.id}.bam
    """
}