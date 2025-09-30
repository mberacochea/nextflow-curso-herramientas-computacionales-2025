process SAMTOOLS_VIEW {
    cpus 2
    memory '4.GB'

    input:
    path sam

    output:
    path "aligned.bam", emit: bam

    script:
    """
    samtools view -@ ${task.cpus} -bS ${sam} > aligned.bam
    """
}