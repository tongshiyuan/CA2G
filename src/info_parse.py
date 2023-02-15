import copy
import pandas as pd
from src.basic_fun import log
from src.parse_fun import parse_cdna, parse_cdna_var


# gene transcript cDNA AA
# MTOR NM_004958 7255G>A E2419K
# IDH1 . 356G>A R119Q
# NM_030649: p.A780T

class var_info:
    def __init__(self, info):
        # chr position ref alt gene transcript cDNA amino_acid_position exon_num
        info_list = info.split('\t')
        self.gene = info_list[0]
        self.transcript = info_list[1]
        self.cdna = info_list[2]
        self.aa = info_list[3].rstrip()
        self.chrom = ''
        self.position = ''
        self.ref = ''
        self.alt = ''
        self.strand = ''
        self.exon_num = ''

    def update(self, tran_info):
        self.chrom = tran_info['chrom']
        self.strand = tran_info['strand']
        self.ref_info = list(tran_info)


def check_var(var, df):
    var_list = []
    if var.transcript and var.transcript != '.':
        info = df[df['transcript'] == var.transcript]
        if info.shape[0] == 1:
            for index, row in info.iterrows():
                var.update(row)
                var_list.append(copy.copy(var))
    elif var.gene and var.gene != '.':
        info = df[df['symbol'] == var.gene]
        for index, row in info.iterrows():
            var.update(row)
            var.transcript = row['transcript']
            var_list.append(copy.copy(var))
    else:
        var_list = 0
    return var_list


def var_parse(file, output, db_file):
    df = pd.read_table(db_file, low_memory=False, header=None)
    df.columns = ['id', 'transcript', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd',
                  'exonCount', 'exonStarts', 'exonEnds', 'score', 'symbol', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
    var_list = []
    with open(file) as f:
        for line in f:
            var = var_info(line)
            rst = check_var(var, df)
            if rst:
                var_list += rst
            else:
                log('WARN', f'Can not found info of <{line}>')
    paesed_var_list = []
    for var in var_list:
        if var.cdna and var.cdna != '.':
            info_list = parse_cdna(var.cdna)
            _var = parse_cdna_var(var, info_list)
            if _var:
                paesed_var_list.append(copy.copy(_var))
        else:
            # 先空着
            pass
    if len(paesed_var_list) > 0:
        fo = open(output, 'w')
        fo.write('chr\tposition\tref\talt\tgene\ttranscript\tstrand\tcDNA\tamino_acid_position\texon_num\n')
        for var in paesed_var_list:
            fo.write('\t'.join([str(i) for i in
                                [var.chrom, var.position, var.ref, var.alt, var.gene, var.transcript, var.strand,
                                 var.cdna, var.aa, var.exon_num]]) + '\n')
        fo.close()
