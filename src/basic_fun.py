import sys
import time


def log(level, text):
    localtime = time.asctime(time.localtime(time.time()))
    if level == 'ERROR':
        sys.exit(f'[{level}]: {localtime} - {text}')
    print(f'[{level}]: {localtime} - {text}')


def read_gene(file='simulator/lib/hg19/refGenePC.txt'):
    transcriptDict = {}
    with open(file) as f:
        for line in f:
            info = line.rstrip().split('\t')
            transcriptDict[info[1]] = {
                'chrom': info[2],
                'strand': info[3],
                'cdsStart': info[6],
                'cdsEnd': info[7],
                'exonStarts': info[9].rstrip(',').split(','),
                'exonEnds': info[10].rstrip(',').split(','),
                'gene': info[12],
                'exonFrames': info[15].rstrip(',').split(','),
            }

    return transcriptDict


def read_seq(file='simulator/lib/hg19/refGeneSeq.txt'):
    seqDict = {}
    with open(file) as f:
        for line in f:
            seqDict[line.split('\t')[0]] = line.split('\t')[1].rstrip()
    return seqDict


def complement(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))


def revcomp(seq):
    return complement(seq)[::-1]
