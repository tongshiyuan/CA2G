import re


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
            _exon = cdna.split('>')[:-3].split('+')[0]
            _intr = cdna.split('>')[:-3].split('+')[1]
            stand = '+'
        elif '-' in cdna:
            _exon = cdna.split('>')[:-3].split('-')[0]
            if not _exon:
                _exon = 0
            _intr = cdna.split('>')[:-3].split('-')[1]
            stand = '-'
        elif '*' in cdna:
            _exon = '*'
            _intr = cdna.split('>')[:-3].split('*')[1]
            stand = '*'
        else:
            _exon = cdna.split('>')[:-3]
            _intr, stand = 0, ''
        return _exon, _intr, stand, ref, alt
    else:
        # 先空着
        pass


def parse_cdna_var(var, info_list):
    # info_list: _exon, _intr, stand, ref, alt
    esl = [int(i) for i in var.ref_info[9].rstrip(',').split(',')]  # exon start list
    eel = [int(i) for i in var.ref_info[10].rstrip(',').split(',')]  # exon end list
    exon_len = [int(j) - int(i) + 1 for i, j in zip(esl, eel)]
    exon_num = len(exon_len)
    tss = int(var.ref_info[6])  # Translation start site
    ted = int(var.ref_info[7])  # Translation end site
    # 判断翻译起始位点的位置
    sl1 = [int(tss >= i) for i in esl]
    sl2 = [int(tss <= i) for i in eel]
    tran_strat_site = [i - j for i, j in zip(sl1, sl2)].index(0)
    el1 = [int(ted >= i) for i in esl]
    el2 = [int(ted <= i) for i in eel]
    tran_end_site = [i - j for i, j in zip(el1, el2)].index(0)
    codon_len = [0] * exon_num
    if tran_end_site > tran_strat_site:
        codon_len[tran_strat_site] = eel[tran_strat_site] - tss + 1
        codon_len[tran_end_site] = ted - esl[tran_end_site] + 1
        for idx in range(tran_strat_site + 1, tran_end_site):
            codon_len[idx] = eel[idx] - esl[idx] + 1
    else:
        codon_len[tran_strat_site] = ted - tss + 1
    std = info_list[2]
    if var.stand == '+':
        if not std:
            for idx, v in enumerate(codon_len):
                tmp_v = info_list[0] - v
                if tmp_v <= 0:
                    pos = eel[idx] + tmp_v
                    break
        else:
            # 先空着
            pass
    else:
        if not std:
            for idx, v in enumerate(codon_len[::-1]):
                tmp_v = info_list[0] - v
                if tmp_v <= 0:
                    pos = esl[exon_num - idx] - tmp_v
                    break
        else:
            # 先空着
            pass
    var.position = pos
    var.exon_num = idx
    var.ref = info_list[3]
    var.alt = info_list[4]
    return var
