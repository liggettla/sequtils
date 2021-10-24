#!/usr/bin/env python3

from multiprocessing import Process, Lock, Queue
import os
import gzip
from Levenshtein import distance
from collections import defaultdict
import pandas as pd

aavs1 = 'GTTAGACCCAATATCAGGAGACTAGGAAGGAGGAGGCCTAAGGATGGGGCTTTTCTGTCACCAATCCTGTCCCTAGTGGCCCCACTGTGGGGTGGAGGGGACAGATAAAAGTACCCAGAACCAGAGCC'
chek2 = 'CTTTATTTCTGCTTAGTGACAGTGCAATTTCAGAATTGTTATTCAAAGGACGGCGTTTTCCTTTCCCTACAAGCTCTGTATTTACAAAGGTTCCATTGCCACTGTGATCTTCTATGTATGCAATGTAAGAGTTTTTAGGACCCACTTCC'
all_edits = defaultdict(lambda: defaultdict(int))

def annotate(filename):
    # filenames and sample numbers
    samples = {
            'TAAGGCGA_GTACGGAT':1,
            'CGTACTAG_CTGTGAAC':2,
            'AGGCAGAA_TCGTAGCA':3,
            'TCCTGAGC_GTGAGTCT':4,

            'GGACTCCT_GCTTGGAT':5,
            'TAGGCATG_TGAATGCG':6,
            'CTCTCTAC_CGTGTCAT':7,
            'CAGAGAGG_CCGTAGTT':8,

            'TAAGGCGA_TAAGCCTC':9,
            'CGTACTAG_ACGAAGTC':10,
            'AGGCAGAA_ACAACGAC':11,
            'TCCTGAGC_GAAGCAGA':12,

            'GGACTCCT_AAGTCCTC':13,
            'TAGGCATG_ATCGAGGA':14,
            'CTCTCTAC_TCCTGCTT':15,
            'CAGAGAGG_GCCTAGAA':16,

            'TAAGGCGA_AGTTACGC':17,
            'CGTACTAG_TGGATCAC':18,
            'AGGCAGAA_AACTGCAC':19,
            'TCCTGAGC_TTACGGTC':20,

            'GGACTCCT_GACGATCA':21,
            'TAGGCATG_CCTTGAAC':22,
            'CTCTCTAC_ACGAATGG':23,
            'CAGAGAGG_GGCATACT':24,
        }

    # requires python 3.9
    tag = filename.removesuffix('_R1.fastq.gz')
    sample = samples[tag]
    edit = ''
    timepoint = ''
    replicate = ''

    if sample in [1,2,3,4]:
        edit = 'AAVS1'
        timepoint = 1
    if sample in [5,6,7,8]:
        edit = 'CHEK2'
        timepoint = 1
    if sample in [9,10,11,12]:
        edit = 'AAVS1'
        timepoint = 2
    if sample in [13,14,15,16]:
        edit = 'CHEK2'
        timepoint = 2
    if sample in [17,18,19,20]:
        edit = 'AAVS1'
        timepoint = 3
    if sample in [21,22,23,24]:
        edit = 'CHEK2'
        timepoint = 3

    if sample in [1,9,17]:
        replicate = 'A'
    if sample in [5,13,21]:
        replicate = 'A'
    if sample in [2,10,18]:
        replicate = 'B'
    if sample in [6,14,22]:
        replicate = 'B'
    if sample in [3,11,19]:
        replicate = 'C'
    if sample in [7,15,23]:
        replicate = 'C'
    if sample in [4,12,20]:
        replicate = 'D'
    if sample in [8,16,24]:
        replicate = 'D'

    return edit, timepoint, replicate

def readFile(indir, filename, q):

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
                        # will edit the gobal dictionary
                        quantifyEdits(filename, edit, test)
                count += 1

            elif count == 3:
                count += 1

            # reset counting
            elif count == 4:
                count = 1

    # put final DataFrame into queue that is
    # shared across processes so it can be
    # combined at the end

    # might need to put the all_edits dictionary in this method instead of making it global
    q.put(pd.DataFrame(all_edits))
    # q.put(5)

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

def quantifyEdits(filename, deletion, edited_sequence):
    # get actual sample numbers and other info
    edit, timepoint, replicate = annotate(filename)

    # stores the sample information
    sample_id = '{}_{}_{}'.format(edit, timepoint, replicate)
    # stores the edited sequence and the number of edits
    deletion_id = '{}_{}'.format(deletion, edited_sequence)

    all_edits[sample_id][deletion_id] += 1

# dedup
def umiCollapse():
    pass

if __name__ == '__main__':
    indir = '../test'

    # just put each DataFrame into the queue and combine at the end
    q = Queue()
    df_combined = pd.DataFrame()

    # get a count for the size of the queue
    queue_count = 0
    # start a process for each input fastq
    for filename in os.listdir(indir):
        if filename.endswith('R1.fastq.gz'):
            queue_count += 1
            p = Process(target=readFile, args=(indir, filename, q))
            p.start()

    # df_combined = df_combined.append(q.get())
    for i in range(queue_count):
        temp = q.get()
        temp = temp.reset_index().rename(columns={'index':'Seq'})
        df_combined = df_combined.append(temp, ignore_index=True)

    df_combined.to_csv('all_edits.txt')
    p.join()
