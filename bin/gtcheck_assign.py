#!/usr/bin/env python3

## pick out the closes matching pair from output by bcftools gtcheck

DEBUG = False

import sys
import statistics

SMVAL = 1.0E-9

def parse_gtcheck(infh):
    gtchkd = {}
    for lin in infh:
        if lin[:2] == 'DC':
            fld = lin.strip().split('\t')
            if DEBUG:
                print(fld)
            if len(fld) > 3 and 'DC' == fld[0].strip():
                donor_query = fld[1].strip()
                donor_ref = fld[2].strip()
                score = float(fld[3])
                tup = (donor_ref, score)
                try:
                    gtchkd[donor_query].append(tup)
                except KeyError:
                    gtchkd[donor_query] = [tup]
    return gtchkd

def find_matches(gtchkd):
    rsltd = {}
    for k in gtchkd:
        gtchkd[k].sort(key = lambda a: a[1])
        if DEBUG:
            print(gtchkd[k])
        donor0, score0 = gtchkd[k][0]
        score1 = 0
        score_n = 0
        n = 0
        d0 = float(0)
        d1 = float(0)
        sd = float(0)
        m = float(0)
        sd = SMVAL
        if len(gtchkd[k]) > 1:
            donor1, score1 = gtchkd[k][1]
            dd1 = score1 - score0
            scores = []
            for i in range(1, len(gtchkd[k])):
                scores.append(gtchkd[k][i][1])
            m = statistics.mean(scores)
            d0 = m - score0
            d1 = m - score1
            sd = statistics.stdev(scores)
            if sd < SMVAL:
                sd = SMVAL
            if d0 < 1.5*sd:
                print("WARNING: ambiguous assignment", donor0, donor1, score0, score1, m, sd, d0, d1)
                sys.stderr.write("ERROR: {:s} assignment ambiguous.\n".format(donor0))
            score_n = scores[-1]
            n = len(scores)
        z0 = d0/sd
        z1 = d1/sd
        rsltd[k] = (donor0, score0, score1, score_n, n, m, sd, z0, z1, dd1/sd)
    return rsltd

if __name__ == '__main__':
    narg = len(sys.argv)
    if narg != 3:
        sys.exit("usage: {:s} <input file> <output file>".format(sys.argv[0]))

    fnin = sys.argv[1]
    fnout = sys.argv[2]

    if fnin == '-':
        infh = sys.stdin
    else:
        infh = open(fnin, 'r')
    oufh = open(fnout, 'w')
    oufh.write("donor_query,donor_gt,score0,score1,score_n,n,mean,sd,z0,z1\n")
    gtchkd = parse_gtcheck(infh)
    if DEBUG:
        print(gtchkd)
    rsltd = find_matches(gtchkd)
    if DEBUG:
        print(rsltd)
    sys.stdout.write("\nAssignments:\n")
    for k in rsltd:
        (donorid, score0, score1, score_n, n, m, sd, z0, z1, zz1) = rsltd[k]
        if z0 < SMVAL:
            donorid = 'UNASSIGNED'
        sys.stdout.write(
            "{:s} -> {:s} score0 = {:.1f}, score1 = {:.1f}, score_last = {:.1f}, n_samples = {:d}, "
            "mean = {:.1f}, sd = {:.1f}, z0 = {:.1f}, z1 = {:.1f}, zz1 = {:.1f}\n"
            .format(k, donorid, score0, score1, score_n, n, m, sd, z0, z1, zz1)
            )
        oufh.write("{:s},{:s},{:.1f},{:.1f},{:.1f},{:d},{:.1f},{:.1f},{:.1f},{:.1f}\n"
            .format(k, donorid, score0, score1, score_n, n, m, sd, z0, z1)
            )
    oufh.close()
    infh.close()

    sys.exit(0)
