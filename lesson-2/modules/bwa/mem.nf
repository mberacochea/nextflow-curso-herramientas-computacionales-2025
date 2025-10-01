process BWA_MEM {

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