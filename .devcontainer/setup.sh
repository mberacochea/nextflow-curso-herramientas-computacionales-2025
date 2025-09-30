#!/usr/bin/env bash

# Customise the terminal command prompt
printf "export PS1='\\[\\e[3;36m\\]\${PWD#/workspaces/} ->\\[\\e[0m\\] '\n" >> $HOME/.bashrc
export PS1='\[\e[3;36m\]${PWD#/workspaces/} ->\[\e[0m\] '

# Initialize micromamba for bash
eval "$(mamba shell hook --shell bash)"
mamba activate

mamba install -c bioconda samtools bwa qualimap

# Update Nextflow
nextflow self-update
nextflow -version

echo ""
echo "=========================================="
echo "Environment setup complete!"
echo "=========================================="
echo ""
echo "Available lessons:"
echo "  - lesson-1: Basic Nextflow pipeline with modules"
echo "  - lesson-2: Centralized configuration"
echo "  - lesson-3: Software environment management (Docker/Conda)"
echo "  - lesson-4: Publishing results"
echo ""
echo "To get started:"
echo "  cd lesson-1"
echo "  cat README.md"
echo ""

cat /usr/local/etc/vscode-dev-containers/first-run-notice.txt