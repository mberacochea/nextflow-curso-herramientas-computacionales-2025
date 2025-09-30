# Nextflow Intro - Herramientas Computacionales 2025

```Introducción a Nextflow y flujos de trabajo reproducibles```

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/mberacochea/nextflow-curso-herramientas-computacionales-2025?quickstart=1&ref=main)

## Learning Objective

Learn how to write and run a Nextflow pipeline starting from a bash script. This hands-on course guides you through converting a bioinformatics alignment workflow into a Nextflow pipeline that is ready to be executed in _almost_ any platform.

## Course Structure

The course is organized into progressive lessons, each building upon the previous one:

### Lesson 1: Basic Nextflow Pipeline

Convert a bash workflow into a Nextflow pipeline, we will start with separating each individual step of the bash script into its own module. We then stitch them together into a workflow, plus a configuration file to the basics.
Learn the fundamentals of processes, channels, and workflow structure while defining resource requirements directly in modules.

Key concepts: Processes, channels, modules, resource directives.

### Lesson 2: Software Environment Management

Add reproducible software environments using Docker containers and Conda. Learn to use profiles for switching between execution modes.

Key concepts: Docker, Conda, profiles, reproducibility

### Lesson 3: Publishing Results

Organize pipeline outputs using `publishDir` directives to create a structured results directory with the pipeline final set of files.

Key concepts: Publishing results, output organization, pattern matching

### Lesson 4: Centralized Configuration

Move resource requirements (cpus, memory, time) from individual modules to a centralized `nextflow.config` file. This lesson includes an exercise to practice centralizing process resources.

Key concepts: Configuration files, process scopes, centralized resource management

### Lesson 6: Processing Multiple Samples

Adapt the pipeline to process multiple samples using a CSV samplesheet input. Introduces the meta map pattern for tracking sample information and organizing outputs by sample ID.

Key concepts: Meta map, CSV input, channel operators, parallel sample processing

## Workflow

The training uses a bioinformatics workflow for read alignment and quality control:

```bash
bwa index reference.fasta

bwa mem reference.fasta test_1.fastq.gz test_1.fastq.gz > aligned.sam

samtools view -bS aligned.sam > aligned.bam

samtools sort aligned.bam -o aligned.sorted.bam

samtools index aligned.sorted.bam

qualimap bamqc -bam aligned.sorted.bam -outdir qualimap_results
```

## Getting Started

### Using GitHub Codespaces

Click the "Code" button above and select "Create codespace on main". The environment will be automatically configured with all necessary tools.

### Using VS Code Dev Containers Locally

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/)
2. Install [VS Code](https://code.visualstudio.com/) with the [Dev Containers extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)
3. Clone this repository
4. Open in VS Code and select "Reopen in Container"

See `.devcontainer/README.md` for detailed testing instructions.

## Running the Lessons

Each lesson directory contains:
- `README.md` - Learning objectives and key changes
- `main.nf` - Main workflow file
- `nextflow.config` - Configuration file
- `modules/` - Process definitions organized by tool

For lesson 1:

```bash
cd lesson-1
nextflow run main.nf \
  --reference /workspaces/nextflow-curso-herramientas-computacionales-2025/assets/genome.fasta \
  --reads '/workspaces/nextflow-curso-herramientas-computacionales-2025/assets/test_{1,2}.fastq.gz' \
  --outdir results
```

For lessons 2-4, specify an execution profile:

```bash
nextflow run main.nf \
  -profile docker \
  --reference /workspaces/nextflow-curso-herramientas-computacionales-2025/assets/genome.fasta \
  --reads '/workspaces/nextflow-curso-herramientas-computacionales-2025/assets/test_{1,2}.fastq.gz' \
  --outdir results
```

For lesson 6 (multiple samples):

```bash
cd lesson-6
nextflow run main.nf \
  -profile docker \
  --reference /workspaces/nextflow-curso-herramientas-computacionales-2025/assets/genome.fasta \
  --samplesheet samplesheet.csv \
  --outdir results
```

## Tools and Technologies

- Nextflow
- Docker
- Conda/Micromamba
- BWA (alignment)
- Samtools (BAM processing)
- Qualimap (quality control)

## Prerequisites

Basic understanding of:
- Command line / bash scripting
- Bioinformatics workflows
- Docker concepts (helpful but not required)

## Course Format

The course has 5 progressive lessons (numbered 1, 2, 3, 4, 6):
- Clear learning objectives
- Comparison with previous lesson
- Working code examples
- README with explanations
- Lessons 1 and 4 include hands-on exercises with solutions

Progress through the lessons sequentially to build a complete understanding of Nextflow pipeline development.
