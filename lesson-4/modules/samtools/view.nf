process SAMTOOLS_VIEW {
    conda "bioconda::samtools=1.17"
    container "biocontainers/samtools:1.17--hd87286a_2"
    publishDir "${params.outdir}/bam", mode: 'copy', pattern: "*.bam"

    input:
    path sam

    output:
    path "aligned.bam", emit: bam

    script:
    """
    samtools view -@ ${task.cpus} -bS ${sam} > aligned.bam
    """
}