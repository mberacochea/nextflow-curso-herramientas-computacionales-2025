process BWA_MEM {
    conda "bioconda::bwa=0.7.17"
    container "biocontainers/bwa:0.7.17--hed695b0_7"
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: "*.sam"

    input:
    tuple path(reference), path(index)
    tuple val(sample_id), path(reads)

    output:
    path "aligned.sam", emit: sam

    script:
    """
    bwa mem -t ${task.cpus} ${reference} ${reads[0]} ${reads[1]} > aligned.sam
    """
}