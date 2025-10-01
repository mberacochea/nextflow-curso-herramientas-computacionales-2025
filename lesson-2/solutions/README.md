# Lesson 2 Solution

## Changes to Make

Add `conda` and `container` directives to all process modules for reproducible software environments.

### BWA Modules (index.nf, mem.nf)

Add after the process declaration:

```groovy
conda "bioconda::bwa=0.7.17"
container "biocontainers/bwa:0.7.17--hed695b0_7"
```

### Samtools Modules (view.nf, sort.nf, index.nf)

Add after the process declaration:

```groovy
conda "bioconda::samtools=1.17"
container "biocontainers/samtools:1.17--hd87286a_2"
```

### Result

The pipeline can now be executed with either Docker or Conda using the `-profile` flag.
