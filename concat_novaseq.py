#!/usr/bin/env python3

import glob
import gzip

# grab files from lane 1
for l1 in glob.iglob('1*.fastq.gz'):

    l1_split = l1.split('.')
    l2 = ('2_HLG52DRXY.2.{}.unmapped.{}.fastq.gz'.format(l1_split[2], l1_split[4]))

    f1 = gzip.open(l1, 'rb')
    f2 = gzip.open(l2, 'rb')
    f3 = gzip.open('{}_R{}.fastq.gz'.format(l1_split[2], l1_split[4]), 'ab')

    # appending the contents of the second file to the first file
    f3.write(f1.read())
    f3.write(f2.read())

    f1.close()
    f2.close()
    f3.close()
