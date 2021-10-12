#!/usr/bin/env python3

from multiprocessing import Process
import os
import gzip

def readFile(indir, filename):
    with gzip.open('{}/{}'.format(indir, filename), 'rb') as f:
        count = 0
        for line in f:
            if count < 10:
                print(filename)
                print(line)
                count += 1

if __name__ == '__main__':
    indir = '../test'
    for filename in os.listdir(indir):
        if filename.endswith('R1.fastq.gz'):
            # readFile(indir, filename)
            p = Process(target=readFile, args=(indir, filename))
            p.start()
    p.join()
