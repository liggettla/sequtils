#!/usr/bin/env python3

import gzip

samples = {
    'a':['ACAAACGG','GTGTCTTA','alex_1'],
    'b':['CACCACAC','TAGGTCTA','alex_2'],
    'c':['ACCCAGCA','TCGATTAG','alex_3'],
    'd':['ACAGTGGT','ACTTAGCA','alex_4'],
    'e':['ATCACGAC','AAGGTTCA','alex_5'],
    'f':['CAGATCCA','AGAGAACA','alex_6'],
    'WT':['CGATGTAT','AGATCTCG','WT'],
    'Mut4':['TTAGGCAT','AGATCTCG','Mut'],
}

# clear previous files
for i in samples:
    filename = samples[i][2]
    f = open(filename + '_R1.fastq','w')
    f.close()

indirs = [
    '/Users/lliggett/Downloads/alex_yong_seq_10.1.2020/all_lanes/Unindexed_Reads-136304169/FASTQ_Generation_2020-10-01_17_05_56Z-322213893/Undetermined_from_200930_NS500410_0123_AHTM5GBGXG_L001-ds.501888d7ad6646099a6cb484453b7541/Undetermined_S0_L001',
    '/Users/lliggett/Downloads/alex_yong_seq_10.1.2020/all_lanes/Unindexed_Reads-136304169/FASTQ_Generation_2020-10-01_17_05_56Z-322213893/Undetermined_from_200930_NS500410_0123_AHTM5GBGXG_L002-ds.fdca5e659e47434aa8fb21c186014e21/Undetermined_S0_L002',
    '/Users/lliggett/Downloads/alex_yong_seq_10.1.2020/all_lanes/Unindexed_Reads-136304169/FASTQ_Generation_2020-10-01_17_05_56Z-322213893/Undetermined_from_200930_NS500410_0123_AHTM5GBGXG_L003-ds.8ded93b68f9148d89dd56432c676543b/Undetermined_S0_L003',
    '/Users/lliggett/Downloads/alex_yong_seq_10.1.2020/all_lanes/Unindexed_Reads-136304169/FASTQ_Generation_2020-10-01_17_05_56Z-322213893/Undetermined_from_200930_NS500410_0123_AHTM5GBGXG_L004-ds.9fe1c4c1a6f14889b4be2c0643e936c1/Undetermined_S0_L004',
]

# concatenate all lanes together
for indir in indirs:
    one = indir + '_R1_001.fastq.gz'
    two = indir + '_R2_001.fastq.gz'

    with gzip.open(one, 'rb') as r1, gzip.open(two, 'rb') as r2:
        linecount = 1
        for line1, line2 in zip(r1, r2):

            if linecount == 1:
                a1 = line1
                a2 = line2
                linecount += 1
            elif linecount == 2:
                b1 = line1
                b2 = line2
                linecount += 1
            elif linecount == 3:
                c1 = line1
                c2 = line2
                linecount += 1
            else:
                d1 = line1
                d2 = line2
                linecount = 1

                line_list = a1.decode('UTF-8').split()
                if len(line_list) == 2:
                    index1 = line_list[1].split(':')[3].split('+')

                line_list = a2.decode('UTF-8').split()
                if len(line_list) == 2:
                    index2 = line_list[1].split(':')[3].split('+')

                # check that indexes match and sort by sample identity into ouptut files
                for i in samples:
                    x = samples[i][0]
                    y = samples[i][1]
                    filename = samples[i][2]

                    if index1 == index2 == [x,y]:
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
