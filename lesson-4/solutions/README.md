# Lesson 4 Solution

## Changes to Make

Add a `process` configuration block to `nextflow.config` to centralize resource management.

### Configuration to Add

Replace the TODO comment in `nextflow.config` with:

```groovy
// Process resource configuration
process {
    // Default resources for all processes
    cpus = 2
    memory = '4.GB'
    time = '2.h'

    // Process-specific resource allocations
    withName: 'BWA_INDEX' {
        cpus = 4
        memory = '8.GB'
    }

    withName: 'BWA_MEM' {
        cpus = 8
        memory = '16.GB'
    }

    withName: 'SAMTOOLS_SORT' {
        cpus = 4
        memory = '8.GB'
    }
}
```

### Result

Resources are now centrally managed in the config file instead of being hardcoded in each process module. The modules use `${task.cpus}` to reference the configured values.
