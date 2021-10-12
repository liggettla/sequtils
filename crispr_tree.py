#!/usr/bin/env python3

from multiprocessing import Process
import os
import gzip
from Levenshtein import distance

def readFile(indir, filename):
    with gzip.open('{}/{}'.format(indir, filename), 'r') as f:
        count = 1
        for line in f:
            # quality
            if count == 1:
                count += 1

            # sequence line
            elif count == 2:
                gene, seq = checkMatch(line.decode("utf-8").strip())
                if gene:
                    print('{} Match in {}:\n{}'.format(gene, filename, seq))
                count += 1

            elif count == 3:
                count += 1

            # reset counting
            elif count == 4:
                count = 1

# check if line contains AAVS1 or CHEK2 sequence
def checkMatch(line):
    aavs1 = 'GTTAGACCCAATATCAGGAGACTAGGAAGGAGGAGGCCTAAGGATGGGGCTTTTCTGTCACCAATCCTGTCCCTAGTGGCCCCACTGTGGGGTGGAGGGGACAGATAAAAGTACCCAGAACCAGAGCC'
    chek2 = 'CTTTATTTCTGCTTAGTGACAGTGCAATTTCAGAATTGTTATTCAAAGGACGGCGTTTTCCTTTCCCTACAAGCTCTGTATTTACAAAGGTTCCATTGCCACTGTGATCTTCTATGTATGCAATGTAAGAGTTTTTAGGACCCACTTCC'

    # if distance(line[10:20], aavs1[0:15]) <= 2:
    if line[10:25] in aavs1:
        return 'AAVS1', aavs1

    #  if distance(line[10:20], chek2[0:15]) <= 2:
    if line[10:25] in chek2:
        return 'CHEK2', chek2

    else:
        return False, False

if __name__ == '__main__':
    indir = '../test'
    for filename in os.listdir(indir):
        if filename.endswith('R1.fastq.gz'):
            # readFile(indir, filename)
            p = Process(target=readFile, args=(indir, filename))
            p.start()
    p.join()
