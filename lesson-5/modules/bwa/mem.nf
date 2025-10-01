process BWA_MEM {
    tag "${meta.id}"
    conda "bioconda::bwa=0.7.17"
    container "biocontainers/bwa:0.7.17--hed695b0_7"
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