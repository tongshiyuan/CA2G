import re
from src.basic_fun import complement

def parse_cdna(cdna):
    '''
    {'replace': 'c.123A>T', 'c.36+1G>T', 'c.*1G>T', 'c.-1G>T',
     'del': 'c.2052delA', 'c.(4071+1_4072-1)_(5154+1_5155-1)del',
     'ins': 'c.5756_5757insAGG', 'c.37+1_37+2insATC',
     'delins': 'c.6775delinsGA',
     'dup': 'c.6_8dupT'
    }
    '''
    cdna = cdna.lstrip('c.')
    if '>' in cdna:
        ref = cdna.split('>')[0][-1]
        alt = cdna.split('>')[1]
        if '+' in cdna:
            _exon = cdna[:-3].split('+')[0]
            _intr = cdna[:-3].split('+')[1]
            stand = '+'
        elif '-' in cdna:
            _exon = cdna[:-3].split('-')[0]
            if not _exon:
                _exon = 0
            _intr = cdna[:-3].split('-')[1]
            stand = '-'
        elif '*' in cdna:
            _exon = '*'
            _intr = cdna[:-3].split('*')[1]
            stand = '*'
        else:
            _exon = cdna[:-3]
            _intr, stand = 0, ''
        return int(_exon), int(_intr), stand, ref, alt
    else:
        # 先空着
        pass


def parse_cdna_var(var, info_list):
    # info_list: _exon, _intr, stand, ref, alt
    esl = [int(i) for i in var.ref_info[9].rstrip(',').split(',')]  # exon start list
    eel = [int(i) for i in var.ref_info[10].rstrip(',').split(',')]  # exon end list
    exon_len = [int(j) - int(i) + 1 for i, j in zip(esl, eel)]
    exon_nums = len(exon_len)
    tss = int(var.ref_info[6])  # Translation start site
    ted = int(var.ref_info[7])  # Translation end site
    # 判断翻译起始位点的位置
    sl1 = [int(tss >= i) for i in esl]
    sl2 = [int(tss <= i) for i in eel]
    tran_strat_site = [i - j for i, j in zip(sl1, sl2)].index(0)
    el1 = [int(ted >= i) for i in esl]
    el2 = [int(ted <= i) for i in eel]
    tran_end_site = [i - j for i, j in zip(el1, el2)].index(0)
    codon_len = [0] * exon_nums
    if tran_end_site > tran_strat_site:
        codon_len[tran_strat_site] = eel[tran_strat_site] - tss
        codon_len[tran_end_site] = ted - esl[tran_end_site]
        for idx in range(tran_strat_site + 1, tran_end_site):
            codon_len[idx] = eel[idx] - esl[idx]
    else:
        codon_len[tran_strat_site] = ted - tss
    if info_list[0] > sum(codon_len):
        return 0
    std = info_list[2]
    if var.strand == '+':
        if not std:
            tmp_v = info_list[0]
            for idx, v in enumerate(codon_len):
                tmp_v = tmp_v - v
                if tmp_v <= 0:
                    if idx == tran_end_site:
                        pos = ted + tmp_v
                    else:
                        pos = eel[idx] + tmp_v
                    var.position = pos
                    ref = info_list[3]
                    alt = info_list[4]
                    exons = idx + 1
                    break
        else:
            # 先空着
            pass
    else:
        if not std:
            tmp_v = info_list[0]
            for idx, v in enumerate(codon_len[::-1]):
                tmp_v = tmp_v - v
                if tmp_v <= 0:
                    if exon_nums - 1 - idx == tran_strat_site:
                        pos = tss - tmp_v + 1
                    else:
                        pos = esl[exon_nums - 1 - idx] - tmp_v + 1
                    # pos = esl[exon_nums - 1 - idx] - tmp_v + 1
                    var.position = pos
                    ref = complement(info_list[3])
                    alt = complement(info_list[4])
                    exons = idx + 1
                    break
        else:
            # 先空着
            pass
    var.exon_num = exons
    var.ref = ref
    var.alt = alt
    return var
