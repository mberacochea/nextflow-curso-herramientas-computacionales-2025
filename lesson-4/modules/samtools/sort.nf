process SAMTOOLS_SORT {
    cpus 2
    memory '1.GB'
    conda "bioconda::samtools=1.17"
    container "quay.io/biocontainers/samtools:1.17--hd87286a_2"
    publishDir "${params.outdir}/bam", mode: 'copy', pattern: "*.sorted.bam"

    input:
    path bam

    output:
    path "aligned.sorted.bam", emit: sorted_bam

    script:
    """
    samtools sort -@ ${task.cpus} ${bam} \\
    -o aligned.sorted.bam
    """
}