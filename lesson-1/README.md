# Lesson 1: Basic Nextflow Pipeline

## Objective

Convert the following bash workflow into a Nextflow pipeline using DSL2 syntax.

## Workflow in Bash

```bash
bwa index reference.fasta
bwa mem reference.fasta reads_R1.fastq.gz reads_R2.fastq.gz > aligned.sam
samtools view -bS aligned.sam > aligned.bam
samtools sort aligned.bam -o aligned.sorted.bam
samtools index aligned.sorted.bam
```

## Nextflow Implementation

The pipeline is organized into individual modules, each representing a step from the bash workflow. Resource requirements (cpus and memory) are defined directly in each module.

## Exercise

The workflow in `main.nf` is incomplete. Steps 4 and 5 are missing. Your task is to complete the workflow by adding:

1. **Step 4**: Sort BAM file using `SAMTOOLS_SORT`
   - Input: `SAMTOOLS_VIEW.out.bam`
   - Output: `sorted_bam`

2. **Step 5**: Index BAM file using `SAMTOOLS_INDEX`
   - Input: `SAMTOOLS_SORT.out.sorted_bam`
   - Output: `indexed_bam`

All modules are already implemented in the `modules/` directory. You only need to add the process calls to the workflow and connect them using the output channels.

### Hints

- Look at how Step 3 uses `SAMTOOLS_VIEW(BWA_MEM.out.sam)`
- Each process output is accessed via `.out.<output_name>`
- The modules for steps 4-5 are in:
  - `modules/samtools/sort.nf`
  - `modules/samtools/index.nf`

## Solution

A complete working solution is available in the `solutions/` directory.

## Running the Pipeline

```bash
nextflow run main.nf \
  --reference /workspaces/nextflow-curso-herramientas-computacionales-2025/assets/genome.fasta \
  --reads '/workspaces/nextflow-curso-herramientas-computacionales-2025/assets/test_{1,2}.fastq.gz' \
  --outdir results
```