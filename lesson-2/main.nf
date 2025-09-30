/*
 * Pipeline parameters
 */
params.reference = null
params.reads = null
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
    if (!params.reads) {
        error "Please provide reads with --reads"
    }

    // Create input channels
    reference_ch = Channel.fromPath(params.reference)
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    // Step 1: Index reference genome
    BWA_INDEX(reference_ch)

    // Step 2: Align reads to reference
    BWA_MEM(
        BWA_INDEX.out.index,
        reads_ch
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