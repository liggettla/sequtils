#!/usr/bin/env python3

#########################################################################
# Purpose of this script is to liftover GWAS sumstats from hg38 to hg19 #
#########################################################################

# worked with 30G of reserved mem
# use gwas_py3 conda environment or informatics3
from pyliftover import LiftOver
import pandas as pd
import sys
import os

indir = os.environ['indir']
summarystats = os.environ['summarystats']
converted = os.environ['converted']

# download chain file
# hg38 to hg19
lo = LiftOver('hg38', 'hg19')

# read in sumstats
sumstats = pd.read_csv('{}/{}'.format(indir, summarystats), sep='\t')

# convert coordinates
chrom = lambda x: lo.convert_coordinate(x.CHR, x.POS)[0][0] if len(lo.convert_coordinate(x.CHR, x.POS)) > 0 else 'chr0'
loc = lambda x: lo.convert_coordinate(x.CHR, x.POS)[0][1] if len(lo.convert_coordinate(x.CHR, x.POS)) > 0 else 0

sumstats['Lifted_Chrom'] = sumstats.apply(chrom, axis='columns')
sumstats['Lifted_Loc'] = sumstats.apply(loc, axis='columns')

# drop unmatched data
sumstats = sumstats[sumstats.Lifted_Loc != 0]

# reformat and eliminate columns unnecessary for ldsc
sumstats.drop(columns=['CHR', 'POS'], inplace=True)
sumstats.rename(columns={'Lifted_Chrom':'CHR',
'Lifted_Loc':'POS',
'Allele1':'A1',
'Allele2':'A2',
'BETA':'Beta',
'SE':'se',
'p.value':'P',},
inplace=True)

# remove chr substring in loc
sumstats.CHR = sumstats.CHR.str.slice(start=3,step=1)

# create UK1000 SNP ID formatted ID
sumstats['SNP'] = sumstats.CHR.astype(str) + ':' + sumstats.POS.astype(str) + '_' + sumstats.A1 + '_' + sumstats.A2

# correct order and only columns needed for the ldsc
sumstats = sumstats[['SNP', 'CHR', 'POS', 'A1', 'A2', 'N', 'Beta', 'se', 'P']]

# export ldsc-ready sumstats
sumstats.to_csv('{}/{}'.format(indir, converted),
sep='\t',
index=False,
compression='gzip',)
