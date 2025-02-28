#!/usr/bin/env bash
INPUT_FILE="$@"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
dt=`date +"%Y_%m_%d_%T"`
cp nextflow.nohup.log ./nextflow.nohup_$dt.log2 || echo 'first time running'
# activate Nextflow conda env

# clean up previous run files
rm -f *.log
rm -f nextflow.nohup.PID.txt 

# start Nextflow in background:
export NXF_OPTS="-Xms5G -Xmx5G"
export NXF_ANSI_LOG=false
export NXF_EXECUTOR_LOG=true
export NXF_LOG_LEVEL=DEBUG
export NXF_CLI_LOG=false
CWD1="$PWD"
parentdir="$(dirname "$CWD1")"
# export RUN_ID="${parentdir##*/}"
export RUN_ID="${PWD##*/}"
mkdir $PWD/work || echo 'exists'
mkdir $PWD/work/tmp || echo 'exists'
echo $RUN_ID | nextflow run $SCRIPT_DIR/../.. -profile sanger -entry JUST_DOUBLETS $INPUT_FILE  --nf_ci_loc $PWD -resume 2>&1 | \
    sed -r "s/\x1b\[[0-9;]*m//g" | \
    ts '[%Y-%m-%d %H:%M:%S]'  | \
    grep -v '^\[[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}\]$' > nextflow.nohup.log &
# get process PID 
sleep 1 && export PID=$(pgrep -f "\\-\\-nf_ci_loc $RUN_DIR")
echo $PID > nextflow.nohup.PID.txt
echo "Nextflow PID is $PID (saved in ./nextflow.nohup.PID.txt)" 
echo kill with \"kill $PID\"
echo "check logs files nextflow.nohup.log and .nextflow.log"
