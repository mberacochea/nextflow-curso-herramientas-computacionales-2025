process BWA_INDEX {
    tag "${reference.simpleName}"
    conda "bioconda::bwa=0.7.17"
    container "biocontainers/bwa:0.7.17--hed695b0_7"
    
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
