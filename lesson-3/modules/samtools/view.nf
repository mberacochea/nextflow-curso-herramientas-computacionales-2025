process SAMTOOLS_VIEW {
    cpus 1
    memory '512.MB'
    conda "bioconda::samtools=1.17"
    container "quay.io/biocontainers/samtools:1.17--hd87286a_2"

    input:
    path sam

    output:
    path "aligned.bam", emit: bam

    script:
    """
    samtools view -@ ${task.cpus} -bS ${sam} > aligned.bam
    """
}