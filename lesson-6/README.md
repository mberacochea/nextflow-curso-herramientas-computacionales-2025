# Lesson 6: Processing Multiple Samples with CSV Input

## Objective

Adapt the pipeline to process multiple samples by reading input from a CSV samplesheet. This lesson introduces the meta map concept, which is essential for tracking sample information throughout the pipeline.

## Workflow in Bash

```bash
bwa index reference.fasta
bwa mem reference.fasta reads_R1.fastq.gz reads_R2.fastq.gz > aligned.sam
samtools view -bS aligned.sam > aligned.bam
samtools sort aligned.bam -o aligned.sorted.bam
samtools index aligned.sorted.bam
qualimap bamqc -bam aligned.sorted.bam -outdir qualimap_results
```

## Key Changes from Lesson 4

### Meta Map

All modules now use a `meta` map to track sample information. The meta map is a Groovy map that contains sample metadata:

```groovy
meta = [id: 'sample1']
```

This approach allows you to:
- Track sample information through the pipeline
- Add additional metadata fields (e.g., patient ID, condition, batch)
- Properly join channels based on sample identity

### Module Changes

Modules now expect `tuple val(meta), path(files)` instead of `tuple val(sample_id), path(files)`:

```groovy
process SAMTOOLS_SORT {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), emit: sorted_bam
    ...
}
```

### Samplesheet Input

The pipeline now reads sample information from a CSV file with the following format:

```csv
sample,fastq_1,fastq_2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

### Samplesheet Parsing

The CSV is parsed inline in the workflow using Nextflow operators:

```groovy
reads_ch = Channel
    .fromPath(params.samplesheet, checkIfExists: true)
    .splitCsv(header: true)
    .map { row ->
        def meta = [id: row.sample]
        def reads = [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]
        return [meta, reads]
    }
```

This creates a channel emitting: `[ [id: 'sample1'], [R1.fastq.gz, R2.fastq.gz] ]`

### Output Organization

Results are now organized by sample ID:

```
results/
├── sample1/
│   ├── alignment/
│   ├── bam/
│   └── qualimap/
├── sample2/
│   ├── alignment/
│   ├── bam/
│   └── qualimap/
└── sample3/
    ├── alignment/
    ├── bam/
    └── qualimap/
```

## Benefits

- Process multiple samples in parallel
- Maintain sample identity throughout the pipeline
- Easily scale from one to hundreds of samples
- Add metadata fields for complex experimental designs
- Properly join results from different processes

## Running the Pipeline

```bash
nextflow run main.nf \
  -profile docker \
  --reference /workspaces/nextflow-curso-herramientas-computacionales-2025/assets/genome.fasta \
  --samplesheet samplesheet.csv \
  --outdir results
```

## Samplesheet Format

Create a `samplesheet.csv` file:

```csv
sample,fastq_1,fastq_2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
sample3,/path/to/sample3_R1.fastq.gz,/path/to/sample3_R2.fastq.gz
```

The pipeline will process all samples in parallel, with Nextflow automatically managing resource allocation and execution.