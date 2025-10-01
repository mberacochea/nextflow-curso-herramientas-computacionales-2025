process SAMTOOLS_SORT {
    cpus 2
    memory '1.GB'

    input:
    path bam

    output:
    path "aligned.sorted.bam", emit: sorted_bam

    script:
    """
    samtools sort -@ ${task.cpus} ${bam} -o aligned.sorted.bam
    """
}