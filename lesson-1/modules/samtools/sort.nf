process SAMTOOLS_SORT {
    cpus 4
    memory '8.GB'

    input:
    path bam

    output:
    path "aligned.sorted.bam", emit: sorted_bam

    script:
    """
    samtools sort -@ ${task.cpus} ${bam} -o aligned.sorted.bam
    """
}