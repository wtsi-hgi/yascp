#!/usr/bin/env bash
CWD1="$PWD"
parentdir="$(dirname "$CWD1")"
INPUT_FILE=$1
export RUN_ID="${PWD##*/}"

# export SINGULARITY_CACHEDIR='/software/hgi/containers/yascp'

export NXF_OPTS="-Xms5G -Xmx5G"
export SINGULARITY_TMPDIR=$PWD/work/tmp
export TEMP=$PWD/work/tmp
export TMP_DIR=$PWD/work/tmp

echo press ENTER to NOT fetch containers, otherwise provide writable path:
read varname

if ["$varname" = ''];
    then
        export NXF_SINGULARITY_CACHEDIR='/software/hgi/containers/yascp'
        export SINGULARITY_DISABLE_CACHE=0
    else
        echo Yascp Will fetch the containers and place them in $varname
        export NXF_SINGULARITY_CACHEDIR=$varname
fi

sample="$RUN_ID.yascp"
echo -e "\nSubmitting yascp (https://github.com/wtsi-hgi/yascp) in JUST_CELLTYPES mode with input file $INPUT_FILE"
bsub -R'select[mem>4000] rusage[mem=4000]' -J yascp_celltypes -n 1 -M 4000 -o yascp_celltypes.o -e yascp_celltypes.e -q normal bash /software/hgi/pipelines/yascp_versions/yascp_v1.4/assets/deploy_scripts/nohup_start_nextflow_lsf_celltypes.sh  $INPUT_FILE
echo "Submitted job can be killed with: bkill -J yascp_celltypes"