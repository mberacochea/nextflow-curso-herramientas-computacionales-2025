# This training container is developed and maintained by Seqera and released under an open-source license (CC BY-NC-ND) for the benefit of the community.

FROM ghcr.io/nextflow-io/training:latest
#
# Creative Commons License:
# Copyright 2025, Seqera. All examples and descriptions are licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.

# Install curl if not present
RUN apt-get update && apt-get install -y curl && rm -rf /var/lib/apt/lists/*

# Install micromamba
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba && \
    mkdir -p /opt/micromamba && cp bin/micromamba /opt/micromamba/
