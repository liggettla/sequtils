#!/usr/bin/env python3

from multiprocessing import Process
import os
import gzip
from Levenshtein import distance
aavs1 = 'GTTAGACCCAATATCAGGAGACTAGGAAGGAGGAGGCCTAAGGATGGGGCTTTTCTGTCACCAATCCTGTCCCTAGTGGCCCCACTGTGGGGTGGAGGGGACAGATAAAAGTACCCAGAACCAGAGCC'
chek2 = 'CTTTATTTCTGCTTAGTGACAGTGCAATTTCAGAATTGTTATTCAAAGGACGGCGTTTTCCTTTCCCTACAAGCTCTGTATTTACAAAGGTTCCATTGCCACTGTGATCTTCTATGTATGCAATGTAAGAGTTTTTAGGACCCACTTCC'

def readFile(indir, filename):
    with gzip.open('{}/{}'.format(indir, filename), 'r') as f:
        count = 1
        for line in f:
            # quality
            if count == 1:
                count += 1

            # sequence line
            elif count == 2:
                line = line.decode("utf-8").strip()
                gene, seq = checkMatch(line)
                if gene:
                    edit, ctrl, test = getEdit(gene, line)
                    if edit:
                        print('Match in {}\nCtrl: {}\nTest: {}\nDistance: {}\n'.format(gene, ctrl, test, edit))
                count += 1

            elif count == 3:
                count += 1

            # reset counting
            elif count == 4:
                count = 1

# check if line contains AAVS1 or CHEK2 sequence
def checkMatch(line):
    # 8bp UMI
    # 14 bp linker

    # if distance(line[10:20], aavs1[0:15]) <= 2:
    if line[25:35] in aavs1:
        return 'AAVS1', aavs1

    #  if distance(line[10:20], chek2[0:15]) <= 2:
    if line[25:35] in chek2:
        return 'CHEK2', chek2

    else:
        return False, False

# locate any CRISPR edits
# maybe convert all characters to the same char
# then distance will just be the difference in length?
def getEdit(gene, line):
    # starts 27bp from center of edits
    starts = {'AAVS1':'GATGGGG',
            'CHEK2':'CTGTATT'
            }
    ends = {'AAVS1':'TGGAGGG',
            'CHEK2':'TGTATGC'
            }
    seqs = {
            'AAVS1':aavs1,
            'CHEK2':chek2
            }

    start = starts[gene]
    end = ends[gene]
    seq = seqs[gene]

    if start in line and end in line:
        ctrl = seq[
                seq.index(start):
                seq.index(end)
                ]

        test = line[
                line.index(start):
                line.index(end)
                ]

        return str(distance(ctrl, test)), ctrl, test
    else:
        return (False, False, False)

# dedup
def umiCollapse():
    pass

if __name__ == '__main__':
    indir = '../test'
    for filename in os.listdir(indir):
        if filename.endswith('R1.fastq.gz'):
            p = Process(target=readFile, args=(indir, filename))
            p.start()
    p.join()
