#!/usr/bin/env python3

import gzip

samples = {
    'AAVS1_DNA':['AAGAAGTTCA','AAVS1_DNA'],
    'AAVS1_RNA':['CGGCGCTGGC','AAVS1_RNA'],
    'RUNX1_DNA':['TCGTCAACTT','RUNX1_DNA'],
    'RUNX1_RNA':['CAACTTGGAT','RUNX1_RNA'],
}

# clear previous files
for i in samples:
    filename = samples[i][1]
    f = open(filename + '_R1.fastq','w')
    f.close()

indirs = [
    '../Undetermined_S0_L001',
]

# concatenate all lanes together
for indir in indirs:
    one = indir + '_R1_001.fastq.gz'
    two = indir + '_R2_001.fastq.gz'
    three = indir + '_I1_001.fastq.gz'

    with gzip.open(one, 'rb') as r1, gzip.open(two, 'rb') as r2, gzip.open(three, 'rb') as r3:
        linecount = 1
        for line1, line2, line3 in zip(r1, r2, r3):

            if linecount == 1:
                a1 = line1
                a2 = line2
                a3 = line3
                linecount += 1
            elif linecount == 2:
                b1 = line1
                b2 = line2
                b3 = line3
                linecount += 1
            elif linecount == 3:
                c1 = line1
                c2 = line2
                c3 = line3
                linecount += 1
            else:
                d1 = line1
                d2 = line2
                d3 = line3
                linecount = 1

                line_list = a1.decode('UTF-8').split('+')
                if len(line_list) == 2:
                    index1 = line_list[1][0:10]

                line_list = a2.decode('UTF-8').split('+')
                if len(line_list) == 2:
                    index2 = line_list[1][0:10]

                line_list = a3.decode('UTF-8').split('+')
                if len(line_list) == 2:
                    index3 = line_list[1][0:10]

                # check that indexes match and sort by sample identity into ouptut files
                for i in samples:
                    x = samples[i][0]
                    filename = samples[i][1]

                    if index1 == index2 == index3 == x:
                        f = open(filename + '_R1.fastq','a')
                        f.write(a1.decode('UTF-8'))
                        f.write(b1.decode('UTF-8'))
                        f.write(c1.decode('UTF-8'))
                        f.write(d1.decode('UTF-8'))
                        f.close()

                        f = open(filename + '_R2.fastq','a')
                        f.write(a2.decode('UTF-8'))
                        f.write(b2.decode('UTF-8'))
                        f.write(c2.decode('UTF-8'))
                        f.write(d2.decode('UTF-8'))
                        f.close()

                        f = open(filename + '_I2.fastq','a')
                        f.write(a3.decode('UTF-8'))
                        f.write(b3.decode('UTF-8'))
                        f.write(c3.decode('UTF-8'))
                        f.write(d3.decode('UTF-8'))
                        f.close()
