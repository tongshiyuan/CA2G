#!/usr/bin/python3
# -*- coding: utf-8 -*-
# author: sytong
# email: tongshiyuan@foxmail.com

# IPO 模型：
# I: gene or transcript，cDNA or amino acid position
# P: if gene and not transcript, get possible transcript
# P: if have amino acid position, get possible cDNA, if have cDNA, check with cDNA
# P: from cNDA to genome position and get possible ref and alt
# O: chr position ref alt gene transcript cDNA amino_acid_position exon_num

import os
import shutil
import argparse
from src.basic_fun import log
from src.info_parse import var_parse

ver = 'v 1.0.0 2023-2-13'


def read_arg():
    parser = argparse.ArgumentParser()
    parser.description = f'cDNA position or Amino acid position to genome position. {ver}'
    parser.add_argument('-i', '--file', required=True, help='variants file.')
    parser.add_argument('-o', '--output', required=False, default='',
                        help='Output of the result. [./<input.file>.pos.txt]')
    parser.add_argument('-r', '--ref', required=False, default='hg38', help='reference version [hg38].')
    # parser.add_argument('-t', '--tmp', required=False, default='.TmpCA2G', help='temp directory [./.TmpCA2G].')
    parser.add_argument('-g', '--gene_db', required=False, default='refGene', help='Gene database version [refGene].')
    args = parser.parse_args()
    print(parser.description)
    return args


def check_par(infile, ref, output):
    if not os.path.isfile(infile):
        log('ERROR', f'Can not find <{infile}>.')
    if os.path.isfile(output):
        log('ERROR', f'File <{output}> existed.')
    if ref not in ['hg19', 'hg38']:
        log('ERROR', f'Can not identify <{ref}>.')


def get_output(arg):
    if not arg.output:
        output = arg.file + '.pos.txt'
    else:
        output = arg.output
    return output


# def get_tmp(tmp_dir):
#     if os.path.isdir(tmp_dir):
#         tmp_dir = tmp_dir.rstrip('/') + 'miaomiaomiao'
#     os.makedirs(tmp_dir)
#     return tmp_dir


def get_db_file(arg, sp):
    db_file_name = f'{arg.ref}_{arg.gene_db}.txt'
    db_file = f'{sp}/DAT/{db_file_name}'
    if not os.path.isfile(db_file):
        log('ERROR', f'Can not find <{db_file_name}>. please check <DAT>.')
    return db_file


def main():
    sp = os.path.split(os.path.realpath(__file__))[0]
    arg = read_arg()
    gene_db = get_db_file(arg, sp)
    infile = arg.file
    output = get_output(arg)
    check_par(infile, gene_db, output)
    # tmp_dir = get_tmp(arg.tmp)
    # try:
    var_parse(infile, output, gene_db)
    # finally:
    #     shutil.rmtree(tmp_dir)


if __name__ == '__main__':
    main()
