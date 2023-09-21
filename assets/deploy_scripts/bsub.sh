#!/usr/bin/env bash
CWD1="$PWD"
parentdir="$(dirname "$CWD1")"
INPUT_FILE=$1
export RUN_ID="${PWD##*/}"
export SINGULARITY_TMPDIR=$PWD/work/tmp
export TEMP=$PWD/work/tmp
export TMP_DIR=$PWD/work/tmp
export SINGULARITY_CACHEDIR='/software/hgi/containers/yascp'
export SINGULARITY_DISABLE_CACHE='False'
sample="$RUN_ID"
echo -e "\n Submitting yascp (https://github.com/wtsi-hgi/yascp) with input file $INPUT_FILE"
bsub -R'select[mem>8000] rusage[mem=8000]' -J $sample -n 1 -M 8000 -o $sample.o -e $sample.e -q long bash /software/hgi/pipelines/yascp/assets/deploy_scripts/nohup_start_nextflow_lsf.sh $INPUT_FILE
echo "Submitted job can be killed with: bkill -J $sample"