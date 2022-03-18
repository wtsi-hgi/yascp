#!/usr/bin/env python3

## pick out the closes matching pair from output by bcftools gtcheck

import sys
import statistics

SMVAL = 1.0E-9

def parse_gtcheck(infh):
    gtchkd = {}
    for lin in infh:
        if lin[:2] == 'DC':
            fld = lin.strip().split('\t')
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
        print(gtchkd[k])
        donor0, score0 = gtchkd[k][0]
        dm = 0
        d1 = 0
        sd = SMVAL
        if len(gtchkd[k]) > 1:
            donor1, score1 = gtchkd[k][1]
            d1 = score1 - score0
            scores = []
            for i in range(1, len(gtchkd[k])):
                scores.append(gtchkd[k][i][1])
            m = statistics.mean(scores)
            dm = m - score0
            sd = statistics.stdev(scores)
            if sd < SMVAL:
                sd = SMVAL
            if dm < 1.5*sd:
                print(scores, m, sd, score0, dm)
                sys.stderr.write("ERROR: assignment possibly ambiguous.\n")

        rsltd[k] = (donor0, dm/sd, d1/sd)
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
    oufh.write("donor_query,donor_gt,z_score,z1_score\n")
    gtchkd = parse_gtcheck(infh)
    print(gtchkd)
    rsltd = find_matches(gtchkd)
    print(rsltd)
    sys.stdout.write("\nAssignments:\n")
    for k in rsltd:
        t = rsltd[k]
        sys.stdout.write("{:s} -> {:s} z = {:f} d1 = {:f}\n".format(k, t[0], t[1], t[2]))
        oufh.write("{:s},{:s},{:.1f},{:.1f}\n".format(k, t[0], t[1], t[2]))
    oufh.close()
    infh.close()

    sys.exit(0)
