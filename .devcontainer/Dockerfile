# This training container is developed and maintained by Seqera and released under an open-source license (CC BY-NC-ND) for the benefit of the community.

FROM ghcr.io/nextflow-io/training:latest
#
# Creative Commons License:
# Copyright 2025, Seqera. All examples and descriptions are licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.

# Install tools
RUN conda install -y -n base -c bioconda samtools bwa qualimap && \
    conda clean -afy && \
    rm -rf /opt/micromamba/pkgs/*
