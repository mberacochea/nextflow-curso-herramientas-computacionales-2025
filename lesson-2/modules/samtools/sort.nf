process SAMTOOLS_SORT {
    conda "bioconda::samtools=1.17"
    container "biocontainers/samtools:1.17--hd87286a_2"

    input:
    path bam

    output:
    path "aligned.sorted.bam", emit: sorted_bam

    script:
    """
    samtools sort -@ ${task.cpus} ${bam} -o aligned.sorted.bam
    """
}