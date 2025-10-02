# Lesson 4: Centralized Configuration

## Objective

Move resource definitions (CPU and memory) from individual process modules to a centralized configuration file. This improves maintainability and makes it easier to adjust resources without modifying the process code.

## Workflow in Bash

```bash
bwa index reference.fasta

bwa mem reference.fasta reads_R1.fastq.gz reads_R2.fastq.gz > aligned.sam

samtools view -bS aligned.sam > aligned.bam

samtools sort aligned.bam -o aligned.sorted.bam

samtools index aligned.sorted.bam
```

## Key Changes from Lesson 3

### Centralized Resource Management

Instead of defining resources in each process module:

```groovy
process BWA_MEM {
    cpus 8
    memory '16.GB'
    ...
}
```

Resources are now managed centrally in `nextflow.config`:

```groovy
process {
    cpus = 2
    memory = '4.GB'
    time = '2.h'

    withName: 'BWA_MEM' {
        cpus = 8
        memory = '16.GB'
    }
}
```

### Exercise

Complete the `nextflow.config` file by adding a `process` configuration block that:

1. Sets default resources for all processes:
   - cpus = 2
   - memory = '4.GB'
   - time = '2.h'

2. Configures process-specific resources:
   - BWA_INDEX: 4 CPUs, 8 GB memory
   - BWA_MEM: 8 CPUs, 16 GB memory
   - SAMTOOLS_SORT: 4 CPUs, 8 GB memory

The modules already reference `${task.cpus}` and will use these configured values automatically.

Check the `solutions/` folder for the complete configuration.

## Benefits

- Single location for resource configuration
- Easy to adjust resources for specific execution environments
- No need to modify process code when tuning resources
- Better separation of concerns (process logic vs. execution configuration)

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
