#!/usr/bin/env python3

## use the donor assgnment file produced by gtcheck_assign_summary.py
## to split the dataset for handover by study

import sys
import os
import csv

FOFN_NAME_INPUT = "files.tsv"
FOFN_NAME_OUTPUT = "files.tsv"
OUTDIR_UNKNOWN = "GT_UNRESOLVED"
STUDY_UNKNOWN = "NONE"

def load_donor_assignments(infn):
    study_dict = {}
    with open(infn, 'r') as infh:
        for row in csv.DictReader(infh):
            study_name = row['panel']
            donor_id = row['donor_query']
            try:
                study_dict[study_name].append(donor_id)
            except KeyError:
                study_dict[study_name] = [donor_id]
    return study_dict

def split_dataset(study_dict, dirnam_input, pool_id, dirnam_output):
    manifest_dict = {}
    # create directories
    for dn in study_dict.keys():
        if dn == STUDY_UNKNOWN:
            d = OUTDIR_UNKNOWN
        else:
            d = dn
        if not os.access(d, os.F_OK):
            os.mkdir(d)
    infh = open(os.path.join(dirnam_input, FOFN_NAME_INPUT), 'r')
    for row in csv.DictReader(infh, delimiter='\t'):
        if pool_id == row['experiment_id']:
            donor_id = row['donor_id']
            for study_name in study_dict:
                if study_name == STUDY_UNKNOWN:
                    dstdn = OUTDIR_UNKNOWN
                else:
                    dstdn = study_name
                dst_dir = os.path.join(dirnam_output, dstdn)
                if donor_id in study_dict[study_name]:
                    src_dir = os.path.join(dirnam_input, pool_id)
                    for fn in (row['filename_h5ad'], row['filename_annotation_tsv']):
                        os.symlink(os.path.join(src_dir, fn), os.path.join(dst_dir, fn))
                    t = (pool_id, donor_id, row['filename_h5ad'], row['filename_annotation_tsv'])
                    try:
                        manifest_dict[dstdn].append(t)
                    except KeyError:
                        manifest_dict[dstdn] = [t]
    infh.close()
    # write manifests
    for dstdn in manifest_dict:
        oufh = open(os.path.join(dirnam_output, dstdn, FOFN_NAME_OUTPUT), 'w')
        oufh.write("pool_id\tdonor_id\tfilename_h5ad\tfilename_annotation_tsv\n")
        for t in manifest_dict[dstdn]:
            oufh.write("{0[0]:s}\t{0[1]:s}\t{0[2]:s}\t{0[3]:s}\n".format(t))
        oufh.close()

    return

if __name__ == '__main__':
    nargs = len(sys.argv)

    if nargs != 4:
        sys.exit(
            "usage: {:s} <donor assignment file [CSV]> <input directory> <pool_id>"
            .format(sys.argv[0])
        )

    fnam_donor_assignments = sys.argv[1]
    dirnam_input = sys.argv[2]
    pool_id = sys.argv[3]
    dirnam_output = os.curdir

    study_dict = load_donor_assignments(fnam_donor_assignments)
    split_dataset(study_dict, dirnam_input, pool_id, dirnam_output)

    sys.exit(0)
