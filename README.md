# Nextflow Intro - Herramientas Computacionales 2025

```Introducción a Nextflow y flujos de trabajo reproducibles```

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/mberacochea/nextflow-curso-herramientas-computacionales-2025?quickstart=1&ref=main)

## Learning Objective

Learn how to write and run a Nextflow pipeline starting from a bash script. This hands-on course guides you through converting a bioinformatics alignment workflow into a Nextflow pipeline that is ready to be executed in _almost_ any platform.

## Course Structure

The course is organized into progressive lessons, each building upon the previous one:

### Lesson 1: Basic Nextflow Pipeline

Convert a bash workflow into a Nextflow pipeline, we will start with separating each individual step of the bash script into its own module. We then stitch them together into a workflow, plus a configuration file for the basics.
Learn the fundamentals of processes, channels, and workflow structure while defining resource requirements directly in modules.

> Key concepts: Processes, channels, modules, resource directives.

### Lesson 2: Software Environment Management

Add reproducible software environments using Docker containers and Conda. Use profiles to switch between execution modes.

> Key concepts: Docker, Conda, profiles, reproducibility

### Lesson 3: Publishing Results

Organize pipeline outputs using publishDir directives to create a structured results directory.

> Key concepts: Publishing results, output organization, pattern matching

### Lesson 4: Centralized Configuration

Move resource requirements (cpus, memory, time) from individual modules to a centralized nextflow.config file.

> Key concepts: Configuration files, process scopes, centralized resource management

### Lesson 5: Processing Multiple Samples

Adapt the pipeline to process multiple samples using a CSV samplesheet input. Introduces the meta map pattern for tracking sample information and organizing outputs by sample ID.

> Key concepts: Meta map, CSV input, channel operators, parallel sample processing

## Workflow

The training uses a bioinformatics workflow for read alignment:

```bash
bwa index reference.fasta

bwa mem reference.fasta test_1.fastq.gz test_1.fastq.gz > aligned.sam

samtools view -bS aligned.sam > aligned.bam

samtools sort aligned.bam -o aligned.sorted.bam

samtools index aligned.sorted.bam
```

Progress through the lessons sequentially to build a complete understanding of Nextflow pipeline development.

---

This training material was developed with assistance from Claude Code.

## Other training resources

The official Nextflow training has loads of high quality courses [https://training.nextflow.io/](https://training.nextflow.io/)