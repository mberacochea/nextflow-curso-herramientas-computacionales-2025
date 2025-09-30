process BWA_INDEX {
    tag "${reference.simpleName}"
    cpus 4
    memory '8.GB'

    input:
    path reference

    output:
    tuple path(reference), path("${reference}.*"), emit: index

    script:
    """
    bwa index -t ${task.cpus} ${reference}
    """
}