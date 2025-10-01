# Lesson 3 Solution

## Changes to Make

Add `publishDir` directives to organize pipeline outputs into structured directories.

### BWA_INDEX (modules/bwa/index.nf)

Already has publishDir - no changes needed.

### BWA_MEM (modules/bwa/mem.nf)

Add after the container directive:

```groovy
publishDir "${params.outdir}/alignment", mode: 'copy', pattern: "*.sam"
```

### SAMTOOLS_VIEW (modules/samtools/view.nf)

Add after the container directive:

```groovy
publishDir "${params.outdir}/bam", mode: 'copy', pattern: "*.bam"
```

### SAMTOOLS_SORT (modules/samtools/sort.nf)

Add after the container directive:

```groovy
publishDir "${params.outdir}/bam", mode: 'copy', pattern: "*.sorted.bam"
```

### SAMTOOLS_INDEX (modules/samtools/index.nf)

Add after the container directive:

```groovy
publishDir "${params.outdir}/bam", mode: 'copy', pattern: "*.bai"
```

### Result

Pipeline outputs will be organized into:
- `${params.outdir}/bwa_index/` - Reference genome index files
- `${params.outdir}/alignment/` - SAM alignment files
- `${params.outdir}/bam/` - BAM files (converted, sorted, and indexed)
