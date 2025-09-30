# Lesson 2: Software Environment Management

## Objective

Add software environment management using Docker containers and Conda to ensure reproducibility across different computing environments.

## Workflow in Bash

```bash
bwa index reference.fasta
bwa mem reference.fasta reads_R1.fastq.gz reads_R2.fastq.gz > aligned.sam
samtools view -bS aligned.sam > aligned.bam
samtools sort aligned.bam -o aligned.sorted.bam
samtools index aligned.sorted.bam
qualimap bamqc -bam aligned.sorted.bam -outdir qualimap_results
```

## Key Changes from Lesson 1

### Container and Conda Directives

Each module now includes both `conda` and `container` directives for software management:

```groovy
process BWA_INDEX {
    conda "bioconda::bwa=0.7.17"
    container "biocontainers/bwa:0.7.17--hed695b0_7"
    ...
}
```

### Execution Profiles

The `nextflow.config` includes profiles for different execution modes:

```groovy
profiles {
    docker {
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
    }

    conda {
        conda.enabled = true
    }
}
```

## Benefits

- Reproducible software environments across different systems
- No need to manually install tools
- Version-controlled software dependencies
- Easy switching between Docker and Conda

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