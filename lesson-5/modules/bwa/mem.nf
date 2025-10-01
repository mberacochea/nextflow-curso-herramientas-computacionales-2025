process BWA_MEM {
    tag "${meta.id}"
    cpus 2
    memory '1.GB'
    conda "bioconda::bwa=0.7.17"
    container "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
    publishDir "${params.outdir}/${meta.id}/alignment", mode: 'copy', pattern: "*.sam"

    input:
    tuple path(reference), path(index)
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.sam"), emit: sam

    script:
    """
    bwa mem -t ${task.cpus} ${reference} ${reads[0]} ${reads[1]} > ${meta.id}.sam
    """
}