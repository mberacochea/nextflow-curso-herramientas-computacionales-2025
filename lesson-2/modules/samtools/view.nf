process SAMTOOLS_VIEW {
    cpus 1
    memory '512.MB'

    input:
    path sam

    output:
    path "aligned.bam", emit: bam

    script:
    """
    samtools view -@ ${task.cpus} -bS ${sam} > aligned.bam
    """
}