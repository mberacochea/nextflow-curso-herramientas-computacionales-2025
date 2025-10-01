process SAMTOOLS_INDEX {
    conda "bioconda::samtools=1.17"
    container "biocontainers/samtools:1.17--hd87286a_2"

    input:
    path sorted_bam

    output:
    tuple path(sorted_bam), path("${sorted_bam}.bai"), emit: indexed_bam

    script:
    """
    samtools index -@ ${task.cpus} ${sorted_bam}
    """
}