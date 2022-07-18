#!/usr/bin/env python3

## make donor assignment across panels

import sys
import csv
import re
import statistics

VERBOSE = True
SMVAL = float(1e-9)
ZSCORE_THRESH = float(8)
ZSCORE_DIST_THRESH = float(8)

DONOR_CELLINE_MATCHSTR = re.compile("^celline_(\S+)$")
FNAM_MATCHSTR = re.compile("^pool_(\S+)_panel_(\S+)_gtcheck_donor_assignments")

class AssignmentTables:
    def __init__(self):
        self.donors = []
        self.cell_lines = {}
        self.panels = {}
        self.cell_line_panel = None

    def parse_assignment(self, infn):
        tabd = {}
        ctr = 0
        fm = FNAM_MATCHSTR.match(infn)
        if not fm:
            sys.exit("ERROR: unexpected input file name '{:s}'.".format(infn))
        panel = fm.group(2)
        infh = open(infn, 'r')
        for row in csv.DictReader(infh):
            ctr += 1
            donor_id = row['donor_query']
            if donor_id not in self.donors:
                self.donors.append(donor_id)
            tabd[row['donor_query']] = (row['donor_gt'], row['score0'], row['score1'], row['z0'], row['z1'])
        infh.close()
        self.panels[panel] = tabd
        return ctr

    def parse_assignment_output_files(self, infns = []):
        ctr = 0
        for fn in infns:
            ctr += 1
            if VERBOSE:
                sys.stderr.write("# input file {:d}: {:s}\n".format(ctr, fn))
            self.parse_assignment(fn)
        spike_in_panel = self.identify_cell_line_panel()
        if VERBOSE:
            sys.stderr.write("# cell line panel is '{:s}'\n".format(spike_in_panel))
        return ctr

    def identify_cell_line_panel(self):
        fnam_cell_line = None
        for nam in self.panels:
            panel = self.panels[nam]
            is_cell_line = None
            is_error = False
            for donor in panel:
                dtup = panel[donor]
                assigned_donor = dtup[0]
                score = float(dtup[1])
                score1 = float(dtup[2])
                cm = DONOR_CELLINE_MATCHSTR.match(assigned_donor)
                if cm:
                    if is_cell_line is None:
                        is_cell_line = True
                        if fnam_cell_line is not None and fnam_cell_line != nam:
                            sys.exit("ERROR: there are multiple panels with cell line names in them, e.g. {:s} and {:s}."
                                .format(fnam_cell_line, fnam))
                        fnam_cell_line = nam
                    else:
                        is_error = not is_cell_line
                    cell_line = cm.group(1)
                    if cell_line in self.cell_lines:
                        self.cell_lines[cell_line][donor] = (score, score1)
                    else:
                        self.cell_lines[cell_line] = {donor: (score, score1)}
                elif is_cell_line is None:
                    is_cell_line = False
                else:
                    is_error = is_cell_line
                if is_error:
                    sys.exit("ERROR: panel {:s} has cell lines mixed with non-cell line labels.".format(panel))
        self.cell_line_panel = fnam_cell_line
        return fnam_cell_line

    def make_cell_line_assignment(self, donor):
        # calculate mean of scores across donors with 1 donor removed
        assignment = None # dictionary of cell lines, values are true/false for a confiden assignment outcome
        n_confident = 0
        for cl in self.cell_lines:
            scores = []
            donor_score = 0
            donor_score1 = 0
            scd = self.cell_lines[cl]
            is_present = False
            for d in self.donors:
                if d not in scd:
                    continue
                sc, sc1 = scd[d]
                if d == donor:
                    donor_score = float(sc)
                    donor_score1 = float(sc1)
                    is_present = True
                else:
                    scores.append(float(sc))

            ns = len(scores)
            if ns > 0:
                m = statistics.mean(scores)
                if ns > 1:
                    sd = statistics.stdev(scores)
                    z = (m - donor_score)/sd
                    is_confident_assignment = z > ZSCORE_THRESH
                    if donor == 'donor8':
                        print("donor8", donor_score, z, sd, m, is_confident_assignment)
                else:
                    is_confident_assignment = 2*donor_score < m
            else:
                is_confident_assignment = donor_score*3 < donor_score1

            if is_confident_assignment and donor_score > SMVAL:
                n_confident += 1
                assignment = cl
        # check that there is only one confident assingment
        if donor == 'donor8':
            print("cell line assignment donor8", donor_score, donor_score1,is_confident_assignment, n_confident)

        if n_confident != 1:
            assignment = None
        return assignment

    def make_panel_assignment(self, panel, donor):
        donor_assigned, score0, score1, z0, z1 = self.panels[panel][donor]
        is_confident_assignment = float(z0) > ZSCORE_THRESH and float(z0)-float(z1) > ZSCORE_DIST_THRESH
        if is_confident_assignment:
            assignment = donor_assigned
        else:
            assignment = None
        return assignment

    def assign_donors(self, oufh):
        oufh.write("donor_query,donor_gt,panel\n")
        if self.cell_line_panel is None:
            self.identify_cell_line_panel()
        for donor in self.donors:
            final_assignment = None
            panel_assignment = None
            cell_line_assignment = None
            panel_ass_panel = None
            cell_line_panel = None
            final_panel = None
            n_assignments = 0
            for panel in self.panels:
                if panel == self.cell_line_panel:
                    cell_line_assignment = self.make_cell_line_assignment(donor)
                    cell_line_panel = panel
                    print("donor:", donor,"panel:", panel, "cell_line_assignment", cell_line_assignment)
                else:
                    assignment = self.make_panel_assignment(panel, donor)
                    print("donor:", donor, ",panel:", panel, "assignment:", assignment)
                    if assignment is not None:
                        panel_assignment = assignment
                        panel_ass_panel = panel
                        n_assignments += 1
            if n_assignments == 1:
                final_assignment = panel_assignment
                final_panel = panel_ass_panel
            elif n_assignments == 0:
                final_assignment = cell_line_assignment
                final_panel = cell_line_panel
            if final_assignment is None:
                final_assignment = 'NONE'
                final_panel = 'NONE'
            if final_panel is None:
                final_panel = 'NONE'
            oustr = "{:s},{:s},{:s}".format(donor, final_assignment, final_panel)
            oufh.write(oustr + '\n')
            print(oustr)
        return

if __name__ == '__main__':
    nargs = len(sys.argv)
    if nargs < 3:
        sys.exit("usage: {:s} <summary output> <assignment output panel 1> [<assignment output panel 2> ....]"
            .format(sys.argv[0]))

    oufn = sys.argv[1]
    infns = sys.argv[2:]

    sys.stdout.write("compare scores from {:d} files and write overall assignment to file {:s}\n".
        format(len(infns), oufn))
    asstab = AssignmentTables()
    fctr = asstab.parse_assignment_output_files(infns)
    if VERBOSE:
        sys.stderr.write("# {:d} assignment files read.\n".format(fctr))
    print(asstab.cell_lines)
    oufh = open(oufn, 'w')
    asstab.assign_donors(oufh)
    oufh.close()

    sys.exit()
