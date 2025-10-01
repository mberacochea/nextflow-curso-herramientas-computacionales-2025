process BWA_INDEX {
    tag "${reference.simpleName}"
    cpus 2
    memory '1.GB'
    conda "bioconda::bwa=0.7.17"
    container "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
    publishDir "${params.outdir}/bwa_index", mode: 'copy'

    input:
    path reference

    output:
    tuple path(reference), path("${reference}.*"), emit: index

    script:
    """
    bwa index -t ${task.cpus} ${reference}
    """
}