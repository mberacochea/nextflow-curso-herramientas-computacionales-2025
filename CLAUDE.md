# Claude Code Project Instructions

## Project Overview

This is a Nextflow training course that teaches how to convert bash bioinformatics workflows into production-ready Nextflow pipelines. The course uses progressive lessons that build on each other.

## Code Style & Conventions

### Nextflow

- Use DSL2 syntax exclusively
- Organize processes into individual module files under `modules/`
- Use simple, hardcoded filenames for single-sample workflows (e.g., `aligned.sam`, `aligned.bam`)
- Only use meta maps in lesson-5 for multiple sample processing
- Always include `${task.cpus}` in tool commands for thread/CPU parameters
- Follow nf-core coding practices

### Configuration Files

**Standard nextflow.config structure:**
```groovy
params {
    reference = null
    reads = null  // or samplesheet for lesson-5
    outdir = "results"
}

// Process resources (lessons 4 and 5 only)
process {
    cpus = 2
    memory = '4.GB'
    time = '2.h'

    withName: 'PROCESS_NAME' {
        cpus = 4
        memory = '8.GB'
    }
}

// Execution reports (lessons 3-5)
timeline {
    enabled = true
    file = "${params.outdir}/timeline.html"
}

// Profiles (lessons 2-5)
profiles {
    docker {
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
    }

    conda {
        conda.enabled = true
        conda.useMicromamba = true
    }
}
```

### Module Structure

**Process modules:**
- Include `conda` and `container` directives (lessons 2+)
- Include `publishDir` for outputs (lessons 3+)
- Use `tag` for process identification when using meta maps
- Keep resource directives in config, not in modules (except lesson-1)

**Example:**
```groovy
process BWA_MEM {
    conda "bioconda::bwa=0.7.17"
    container "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: "*.sam"

    input:
    tuple path(reference), path(index)
    tuple val(sample_id), path(reads)

    output:
    path "aligned.sam", emit: sam

    script:
    """
    bwa mem -t ${task.cpus} ${reference} ${reads[0]} ${reads[1]} > aligned.sam
    """
}
```

## Documentation Style

### README Files

- Use minimal emphasis - avoid excessive bold formatting
- Start objectives with action verbs (e.g., "Add software environment management...")
- Keep sections concise and direct
- Include bare bash workflow at the top for context
- Always show both Docker and Conda running examples (lessons 2+)

### Exercise Lessons

**Lessons with exercises (1, 2, 3, 4):**
- Create `solutions/README.md` with instructions on what changes to make
- Do NOT copy entire code files to solutions (except lesson-1 and lesson-4 config)
- Use clear, step-by-step instructions in solutions/README.md

## Testing & Running

**Standard test data paths:**
```bash
--reference /workspaces/nextflow-curso-herramientas-computacionales-2025/assets/genome.fasta
--reads '/workspaces/nextflow-curso-herramientas-computacionales-2025/assets/test_{1,2}.fastq.gz'
--outdir results
```

## Tools & Environment

- Use conda for dependency management (not mamba install directly in Dockerfile)
- Container: Base on `ghcr.io/nextflow-io/training:latest`
- Tools: samtools, bwa (no qualimap - removed to keep container lean)
- GitHub Codespaces enabled with devcontainer setup

## Important Notes

- No qualimap in this course (removed for container size)
- Workflow ends at samtools index (5 steps total)
- Only lesson-5 processes multiple samples
- Keep configurations consistent across all lessons
- Use profiles for execution, never hardcode docker.enabled outside profiles
