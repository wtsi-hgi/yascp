#!/usr/bin/env python
import argparse
import sys
from pathlib import Path
import gzip
from scipy.stats.stats import pearsonr
import pandas
import seaborn
import matplotlib.pyplot
import numpy
import glob

if sys.version_info < (3, 5):
    print("Python version:")
    print(sys.version_info)
    raise RuntimeError("Requires Python 3.x")

parser = argparse.ArgumentParser(description='Match donors from multiple Vireo VCF outputs.')
parser.add_argument('cellSNP_dirs',
                    help="Directory containing cellSNP output dirs. Subdirectories in this directory " +
                         "should contain GT_donors.vireo.vcf.gz. First directory name is taken as sample name.")
parser.add_argument('outdir', help='directory to write the matcher outfiles')
parser.add_argument('-m','--min_cor_threshold', help='The minimum correlation at which two donors are considered the same', default=0.8)
parser.add_argument('-n','--nonverbose', help='Do not print all info to console',
                    action='store_true')
parser.add_argument('-s','--samplesheet', help='file with mapping between IDs. If this is used, ID is the 4th column')

args = parser.parse_args()
sampleid_to_name = {}

print('Search for Vireo files in subdirectories of: '+args.cellSNP_dirs)
if args.samplesheet:
    print('Use 1st and 4th column of '+args.samplesheet+' for matching vireo names to sample names')


if args.samplesheet:
    # Matching sample IDs between sanger. Hardcoded first and 4th column because that's where the IDs
    # are in my samplesheet, have to think of how to adjust that
    with open(args.samplesheet) as input_file:
        for line in input_file:
            line = line.strip().split('\t')
            sampleid_to_name[line[0]] = line[3]


donor_genotypes = {}
genotypes_per_file = {}
all_chr_pos = set([])
donors = set([])
# Recursively search the args.cellSNP_dirs directory for files that are called GT_donors.vireo.vcf.gz
# This is where Vireo put the genotypes
x = 0 
for f in glob.glob(f"{args.cellSNP_dirs}/*/*T_donors.vire*.vcf.gz"):
    x += 1
    print(f)
    # if x>5:
    #     break
    if not args.nonverbose:
        print('Reading '+str(f)+'...')
    try:
        # Try to open as a gzipped file
        input_file = gzip.open(f, 'rt')
        line = input_file.readline().strip().split('\t')
    except gzip.BadGzipFile:
        # If not gzipped, open as a regular text file
        input_file = open(f, 'rt')
        line = input_file.readline().strip().split('\t')

    with input_file:
        cellranger_id = str(f).split('/')[-2].replace('vireo_','').replace('gt_rep_','')
        # if cellranger_id in sampleid_to_name:
        #     name = sampleid_to_name[cellranger_id]
        # else:
        name = cellranger_id
        # Get the name I use from the cellranger ID through the samplesheet provided
        for line in input_file:
            line = line.strip().split('\t')
            if line[0].startswith('##'):
                continue
            if line[0].startswith('#CHROM'):
                donor_IDs = line[9:]
                # File is a VCF format so donor IDs can be retrieved from #CHROM line
                for donor in donor_IDs:
                    if (donor=='doublet'):
                        continue
                    # The donor is added to the sample so each sample.donor combination will have a list of genotypes
                    donor_genotypes[name+'<SPLIT>'+donor] = []
                    donors.add(name+'<SPLIT>'+donor)
                continue
            gt_index = line[8].split(':').index('GT')
            chr = line[0]
            pos = line[1]
            if chr+'_'+pos not in genotypes_per_file:
                # Use the chr+pos to save the genotype
                genotypes_per_file[chr+'_'+pos] = {}

            all_chr_pos.add(chr+'_'+pos)
            for index, donor in enumerate(line[9:]):
                alleles = donor.split(':')[gt_index]
                alleles = alleles.replace('|',"/")
                # print(alleles)
                # for later comparison reasons, encode the alleles as dosages
                if alleles == '1/0':
                    genotype = 1
                elif alleles == '0/1':
                    genotype = 1
                elif alleles == '0/0':
                    genotype = 0
                elif alleles == '1/1':
                    genotype = 2
                elif alleles == './.':
                    continue
                else:
                    print('Unknown allele configuration: '+alleles)
                    continue
                    # raise RuntimeError('Unknown allele configuration: '+alleles)
                genotypes_per_file[chr + '_' + pos][name+'<SPLIT>'+donor_IDs[index]] = genotype

print('Read',x,'genotype files')
if x == 0:
    print('ERROR: Becasue',x,'files where read, exiting now')
    exit()

# loop over all the observed chr and positions
for chr_pos in all_chr_pos:
    # loop over all the donors
    for donor in donor_genotypes:
        # If this chr.pos was not seen for this specific sample.donor, add None to the list
        # Otherwise, add the genotype. This should make equal lists for all samples with
        # missing values if genotype was not observed
        # These can later be compared
        if donor not in genotypes_per_file[chr_pos]:
            donor_genotypes[donor].append(None)
        else:
            donor_genotypes[donor].append(genotypes_per_file[chr_pos][donor])


donors = list(donors)
correlation_matrix = []
donors_matching = {}
max_group = 0
for donor1 in donors:
    print('Correlating donor '+donor1)
    cor_list = []
    for donor2 in donors:
        tmp_don1 = []
        tmp_don2 = []
        # We want to use only those genotypes that have been called by both donor1 and donor2 (so not none in either of the lists)
        # Make a new, temporary list with only genotypes called by both
        for index in range(0, len(donor_genotypes[donor1])):
            if donor_genotypes[donor1][index] != None and donor_genotypes[donor2][index] != None:
                tmp_don1.append(donor_genotypes[donor1][index])
                tmp_don2.append(donor_genotypes[donor2][index])
        # Calculate correlation between the two lists where all values where one or the other was None have been removed
        # Because we reannotated genotypes to dosages we can correlate directly, so that A/T T/T is better than A/A T/T (because A/A=0, A/T=1, T/T = 2)
        # Think this works better than exact match because with RNA if mistake is made it's most often het -> ref hom
        pcor = pearsonr(tmp_don1,tmp_don2)[0]
        if pcor > float(args.min_cor_threshold):
            # Have to add donors to the same group if one of the donors was already added earlier in loop, so check all existing groups
            # if one of the donors is in it, if so, save which group this is
            donor_already_in_group = False

            for group in donors_matching:
                if donor1 in donors_matching[group] or donor2 in donors_matching[group]:
                    donor_already_in_group = True
                    group_counter = group

            if not donor_already_in_group:
                group_counter = max_group+1

            if group_counter not in donors_matching:
                donors_matching[group_counter] = set([])
            if not args.nonverbose:
                print('Correlation > '+str(args.min_cor_threshold)+' for '+donor1+' and '+donor2+'. Add to group '+str(group_counter))
            donors_matching[group_counter].add(donor1)
            donors_matching[group_counter].add(donor2)

            if max_group < group_counter:
                max_group = group_counter

        cor_list.append(pcor)
    correlation_matrix.append(cor_list)

group_count = 0
with open(args.outdir+'/matched_donors.txt','w') as out:
    out.write('new_donor\tsample\told_donor\n')
    for group in donors_matching:
        group_count += 1
        for sample in donors_matching[group]:
            out.write('donor'+str(group_count)+'\t'+sample.split('<SPLIT>')[0]+'\t'+sample.split('<SPLIT>')[1]+'\n')

correlation_dataframe = pandas.DataFrame(correlation_matrix)
correlation_dataframe.columns = [x.replace('<SPLIT>',' - ') for x in donors]
correlation_dataframe.index = [x.replace('<SPLIT>',' - ') for x in donors]

# plotting correlation heatmap
import numpy as np
correlation_dataframe = correlation_dataframe.fillna(0)  # Replace NaN with 0
correlation_dataframe = correlation_dataframe.replace([np.inf, -np.inf], 0)  # Replace inf/-inf with 0

dataplot = seaborn.clustermap(correlation_dataframe, cmap="YlGnBu", yticklabels=1, xticklabels=1)
dataplot.ax_heatmap.set_xticklabels(dataplot.ax_heatmap.get_xmajorticklabels(), fontsize = 6)
dataplot.ax_heatmap.set_yticklabels(dataplot.ax_heatmap.get_ymajorticklabels(), fontsize = 6)
#mask = numpy.tril(numpy.ones_like(correlation_dataframe))
#values = dataplot.ax_heatmap.collections[0].get_array().reshape(correlation_dataframe.shape)
#new_values = numpy.ma.array(values, mask=mask)
#dataplot.ax_heatmap.collections[0].set_array(new_values)
correlation_dataframe.to_csv('donor_corelations_matrix.tsv',sep='\t')

# displaying heatmap
#matplotlib.pyplot.show()

if not Path(args.outdir).is_dir():
    print(args.outdir+' does not exist, making directory')
Path(args.outdir).mkdir(parents=True, exist_ok=True)
matplotlib.pyplot.savefig(args.outdir+'/correlations.png',dpi=300)

