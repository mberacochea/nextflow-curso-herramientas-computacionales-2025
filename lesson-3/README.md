# Lesson 3: Publishing Results

## Objective

Add `publishDir` directives to modules to organize and save final outputs in a structured results directory.

## Workflow in Bash

```bash
bwa index reference.fasta
bwa mem reference.fasta reads_R1.fastq.gz reads_R2.fastq.gz > aligned.sam
samtools view -bS aligned.sam > aligned.bam
samtools sort aligned.bam -o aligned.sorted.bam
samtools index aligned.sorted.bam
```

## Key Changes from Lesson 2

### publishDir Directives

Each module now includes a `publishDir` directive to save results to organized output directories:

```groovy
process SAMTOOLS_SORT {
    publishDir "${params.outdir}/bam", mode: 'copy', pattern: "*.sorted.bam"
    ...
}
```

### Output Structure

The pipeline now organizes results into subdirectories:

- `${params.outdir}/bwa_index/` - Reference genome index files
- `${params.outdir}/alignment/` - SAM alignment files
- `${params.outdir}/bam/` - BAM files (converted, sorted, and indexed)

### Pattern Matching

The `pattern` parameter in `publishDir` allows selective publishing of specific file types, preventing intermediate files from cluttering the output directory.

## Benefits

- Organized output structure
- Easy to find final results
- Separation of intermediate and final files
- Reproducible directory structure across runs

## Running the Pipeline

With Docker:
```bash
nextflow run main.nf \
  -profile docker \
  --reference /workspaces/nextflow-curso-herramientas-computacionales-2025/assets/genome.fasta \
  --reads '/workspaces/nextflow-curso-herramientas-computacionales-2025/assets/test_{1,2}.fastq.gz' \
  --outdir results
```

With Conda:
```bash
nextflow run main.nf \
  -profile conda \
  --reference /workspaces/nextflow-curso-herramientas-computacionales-2025/assets/genome.fasta \
  --reads '/workspaces/nextflow-curso-herramientas-computacionales-2025/assets/test_{1,2}.fastq.gz' \
  --outdir results
```