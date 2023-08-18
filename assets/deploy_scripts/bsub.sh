#!/usr/bin/env bash
CWD1="$PWD"
parentdir="$(dirname "$CWD1")"
INPUT_FILE=$1
export RUN_ID="${PWD##*/}"
sample="$RUN_ID.yascp"
echo "submitting yascp with input file $INPUT_FILE"
bsub -R'select[mem>4000] rusage[mem=4000]' -G hgi -J $sample -n 1 -M 4000 -o $sample.o -e $sample.e -q long bash /software/hgi/pipelines/yascp/assets/deploy_scripts/nohup_start_nextflow_lsf.sh $INPUT_FILE