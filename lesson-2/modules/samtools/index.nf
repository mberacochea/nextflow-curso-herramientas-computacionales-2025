process SAMTOOLS_INDEX {
    cpus 1
    memory '512.MB'
    conda "bioconda::samtools=1.17"
    container "quay.io/biocontainers/samtools:1.17--hd87286a_2"

    input:
    path sorted_bam

    output:
    tuple path(sorted_bam), path("${sorted_bam}.bai"), emit: indexed_bam

    script:
    """
    samtools index -@ ${task.cpus} ${sorted_bam}
    """
}