import sys
import time


def log(level, text):
    localtime = time.asctime(time.localtime(time.time()))
    if level == 'ERROR':
        sys.exit(f'[{level}]: {localtime} - {text}')
    print(f'[{level}]: {localtime} - {text}')


def complement(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))


def revcomp(seq):
    return complement(seq)[::-1]
