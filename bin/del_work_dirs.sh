#!/usr/bin/env bash
workdir=$1
echo running del_work_dirs_failed.sh script...
rm -r $workdir
rm nextflow.nohup_2*
rm .nextflow.log*