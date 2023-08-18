#!/usr/bin/env bash
CWD1="$PWD"
parentdir="$(dirname "$CWD1")"
INPUT_FILE=$1
export RUN_ID="${PWD##*/}"
sample="$RUN_ID.yascp"
echo "Submitting yascp (https://github.com/wtsi-hgi/yascp) in test mode withsample OneK1k dataset"
bsub -R'select[mem>4000] rusage[mem=4000]' -J yascp_test -n 1 -M 4000 -o yascp_test.o -e yascp_test.e -q long bash /software/hgi/pipelines/yascp/assets/deploy_scripts/nohup_start_nextflow_lsf_test.sh
echo "Submitted job can be killed with: bkill -J yascp_test"