process SAMTOOLS_INDEX {
    cpus 1
    memory '512.MB'

    input:
    path sorted_bam

    output:
    tuple path(sorted_bam), path("${sorted_bam}.bai"), emit: indexed_bam

    script:
    """
    samtools index -@ ${task.cpus} ${sorted_bam}
    """
}