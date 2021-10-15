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
                    print('{} Match in {}:\n{}'.format(gene, filename, line))
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
def getEdit():
    pass
    sub = 'GGGCT'
    ctrl_start = aavs1.index(sub)
    ctrl = (ctrl_start,ctrl_start+20)
    test_start = line.index(sub)
    test = (test_start,test_start+20)
    print(distance(ctrl, test))

if __name__ == '__main__':
    indir = '../test'
    for filename in os.listdir(indir):
        if filename.endswith('R1.fastq.gz'):
            p = Process(target=readFile, args=(indir, filename))
            p.start()
    p.join()
