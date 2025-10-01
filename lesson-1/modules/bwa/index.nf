process BWA_INDEX {
    tag "${reference.simpleName}"
    cpus 2
    memory '1.GB'

    input:
    path reference

    output:
    tuple path(reference), path("${reference}.*"), emit: index

    script:
    """
    bwa index ${reference}
    """
}