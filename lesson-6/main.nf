/*
 * Pipeline parameters
 */
params.samplesheet = null
params.reference = null
params.outdir = "results"

/*
 * Import modules
 */
include { BWA_INDEX      } from './modules/bwa/index'
include { BWA_MEM        } from './modules/bwa/mem'
include { SAMTOOLS_VIEW  } from './modules/samtools/view'
include { SAMTOOLS_SORT  } from './modules/samtools/sort'
include { SAMTOOLS_INDEX } from './modules/samtools/index'
include { QUALIMAP_BAMQC } from './modules/qualimap/bamqc'

/*
 * Main workflow
 */
workflow {
    // Input validation
    if (!params.reference) {
        error "Please provide a reference genome with --reference"
    }
    if (!params.samplesheet) {
        error "Please provide a samplesheet with --samplesheet"
    }

    // Create input channels
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true)

    // Parse samplesheet and create channel with [meta, [reads]]
    reads_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            def reads = [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]
            return [meta, reads]
        }

    // Step 1: Index reference genome
    BWA_INDEX(reference_ch)

    // Step 2: Align reads to reference
    BWA_MEM(
        BWA_INDEX.out.index,
        reads_ch.first()
    )

    // Step 3: Convert SAM to BAM
    SAMTOOLS_VIEW(BWA_MEM.out.sam)

    // Step 4: Sort BAM file
    SAMTOOLS_SORT(SAMTOOLS_VIEW.out.bam)

    // Step 5: Index BAM file
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.sorted_bam)

    // Step 6: Run quality control
    QUALIMAP_BAMQC(SAMTOOLS_INDEX.out.indexed_bam)
}