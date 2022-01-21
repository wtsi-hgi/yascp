#!/usr/bin/env bash

eval "$(conda shell.bash hook)"
conda activate nextflow

export SINGULARITY_CACHEDIR="$PWD/singularity_cache"
mkdir -p $SINGULARITY_CACHEDIR

nextflow run main.nf \
	 -profile sanger \
	 --singularity_pull_docker_container \
	 -c sample_input/gn5_dev/input_run_cellbender.nf \
	 -resume
