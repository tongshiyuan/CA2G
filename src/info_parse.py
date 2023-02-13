from src.basic_info import codon, aaDict
from src.basic_fun import revcomp


def var_parse(file, ref, output, sp):
    pass


def parse_snv(nmid, aa, gene, seqDict, transcriptDict):
    #  NM_030649:p.A780T
    refAA = aa[0]  # A
    altAA = aa[-1]  # T
    pcPos = eval(aa[1:-1])  # 780

    strand = transcriptDict[nmid]['strand']  # + / -
    chrom = transcriptDict[nmid]['chrom']  # chr1
    alt1, alt2, alt3 = [], [], []
    if strand == '+':
        pst = pcPos * 3 - 3  # 0
        ped = pcPos * 3  # 3
        if codon[str.upper(seqDict[nmid][pst:ped])] == refAA:
            ref1 = str.upper(seqDict[nmid][pcPos * 3 - 3])
            ref2 = str.upper(seqDict[nmid][pcPos * 3 - 2])
            ref3 = str.upper(seqDict[nmid][pcPos * 3 - 1])
        else:
            print('<', gene, nmid, aa, '> have wrong reference !')
            return []
        a1, a2, a3 = ['A', 'T', 'C', 'G'], ['A', 'T', 'C', 'G'], ['A', 'T', 'C', 'G']
        a1.remove(ref1)
        for r1 in a1:
            if codon[r1 + ref2 + ref3] == altAA:
                alt1.append(r1)
        a2.remove(ref2)
        for r2 in a2:
            if codon[ref1 + r2 + ref3] == altAA:
                alt2.append(r2)
        a3.remove(ref3)
        for r3 in a3:
            if codon[ref1 + ref2 + r3] == altAA:
                alt3.append(r3)

        cdsStart = eval(transcriptDict[nmid]['cdsStart'])
        cdsEns = eval(transcriptDict[nmid]['cdsEnd'])
        exonStarts = transcriptDict[nmid]['exonStarts']
        exonEnds = transcriptDict[nmid]['exonEnds']
        exonFrames = transcriptDict[nmid]['exonFrames']
        for i in range(len(exonStarts)):
            if eval(exonStarts[i]) > cdsStart:
                sdx = i - 1
                break
        else:
            sdx = len(exonStarts) - 1
        for i in range(len(exonEnds)):
            if eval(exonEnds[i]) >= cdsEns:
                edx = i
                break
        for idx in range(sdx, edx + 1):
            if idx == sdx:
                if pcPos - (eval(exonEnds[idx]) - cdsStart) // 3 > 1:
                    pcPos -= (eval(exonEnds[idx]) - cdsStart) // 3
                elif pcPos - (eval(exonEnds[idx]) - cdsStart) // 3 == 1 and (eval(exonEnds[idx]) - cdsStart) % 3 != 0:
                    if (eval(exonEnds[idx]) - cdsStart) % 3 == 1:
                        gpos1 = eval(exonEnds[idx])
                        gpos2 = eval(exonStarts[idx + 1]) + 1
                        gpos3 = eval(exonStarts[idx + 1]) + 2
                        break
                    elif (eval(exonEnds[idx]) - cdsStart) % 3 == 2:
                        gpos1 = eval(exonEnds[idx]) - 1
                        gpos2 = eval(exonEnds[idx])
                        gpos3 = eval(exonStarts[idx + 1]) + 1
                        break
                else:
                    gpos1 = cdsStart + pcPos * 3 - 2
                    gpos2 = cdsStart + pcPos * 3 - 1
                    gpos3 = cdsStart + pcPos * 3
                    break
            elif idx < edx:
                if (pcPos - (eval(exonEnds[idx]) - eval(exonStarts[idx]) + eval(exonFrames[idx])) // 3) > 1:
                    pcPos -= (eval(exonEnds[idx]) - eval(exonStarts[idx]) + eval(exonFrames[idx])) // 3
                elif pcPos - (
                        eval(exonEnds[idx]) - eval(exonStarts[idx]) + eval(exonFrames[idx])) // 3 == 1 and (
                        eval(exonEnds[idx]) - eval(exonStarts[idx]) + eval(exonFrames[idx])) % 3 != 0:
                    if (eval(exonEnds[idx]) - eval(exonStarts[idx]) + eval(exonFrames[idx])) % 3 == 1:
                        gpos1 = eval(exonEnds[idx])
                        gpos2 = eval(exonStarts[idx + 1]) + 1
                        gpos3 = eval(exonStarts[idx + 1]) + 2
                        break
                    elif (eval(exonEnds[idx]) - eval(exonStarts[idx]) + eval(exonFrames[idx])) % 3 == 2:
                        gpos1 = eval(exonEnds[idx]) - 1
                        gpos2 = eval(exonEnds[idx])
                        gpos3 = eval(exonStarts[idx + 1]) + 1
                        break
                else:
                    gpos1 = eval(exonStarts[idx]) - eval(exonFrames[idx]) + pcPos * 3 - 2
                    gpos2 = eval(exonStarts[idx]) - eval(exonFrames[idx]) + pcPos * 3 - 1
                    gpos3 = eval(exonStarts[idx]) - eval(exonFrames[idx]) + pcPos * 3
                    break
            elif idx == edx:
                if (pcPos - (cdsEns - eval(exonStarts[idx]) + eval(exonFrames[idx])) // 3) > 1 or (pcPos - (
                        cdsEns - eval(exonStarts[idx]) + eval(exonFrames[idx])) // 3 == 1 and (cdsEns - eval(
                    exonStarts[idx]) + eval(exonFrames[idx])) % 3 != 0):
                    print(gene, '超出转录本长度')
                    gpos1, gpos2, gpos3 = 0, 0, 0
                else:
                    gpos1 = eval(exonStarts[idx]) - eval(exonFrames[idx]) + pcPos * 3 - 2
                    gpos2 = eval(exonStarts[idx]) - eval(exonFrames[idx]) + pcPos * 3 - 1
                    gpos3 = eval(exonStarts[idx]) - eval(exonFrames[idx]) + pcPos * 3
                    break

    else:
        pst = -pcPos * 3 - 1  # -4
        ped = -pcPos * 3 + 2  # -1
        if codon[complement(seqDict[nmid][ped:pst:-1])] == refAA:
            ref1 = seqDict[nmid][-pcPos * 3]  # -3
            ref2 = seqDict[nmid][-pcPos * 3 + 1]  # -2
            ref3 = seqDict[nmid][-pcPos * 3 + 2]  # -1
        else:
            print('<', gene, nmid, aa, '> have wrong reference !')
            return []
        a1, a2, a3 = ['A', 'T', 'C', 'G'], ['A', 'T', 'C', 'G'], ['A', 'T', 'C', 'G']
        a1.remove(ref1)
        for r in a1:
            if codon[revcomp(r + ref2 + ref3)] == altAA:
                alt1.append(r)
        a2.remove(ref2)
        for r in a2:
            if codon[revcomp(ref1 + r + ref3)] == altAA:
                alt2.append(r)
        a3.remove(ref3)
        for r in a3:
            if codon[revcomp(ref1 + ref2 + r)] == altAA:
                alt3.append(r)
        cdsStart = eval(transcriptDict[nmid]['cdsEnd'])
        cdsEns = eval(transcriptDict[nmid]['cdsStart'])
        exonStarts = transcriptDict[nmid]['exonEnds']
        exonEnds = transcriptDict[nmid]['exonStarts']
        exonFrames = transcriptDict[nmid]['exonFrames']
        for i in range(len(exonEnds)):
            if eval(exonEnds[i]) > cdsEns:
                edx = i - 1
                break
        else:
            edx = len(exonEnds) - 1

        for i in range(len(exonStarts)):
            if eval(exonStarts[i]) >= cdsStart:
                sdx = i
                break
        for idx in range(sdx, edx - 1, -1):
            if idx == sdx:
                if pcPos - (cdsStart - eval(exonEnds[idx])) // 3 > 1:
                    pcPos -= (cdsStart - eval(exonEnds[idx])) // 3
                elif pcPos - (cdsStart - eval(exonEnds[idx])) // 3 == 1 and (
                        cdsStart - eval(exonEnds[idx])) % 3 != 0:
                    if (cdsStart - eval(exonEnds[idx])) % 3 == 1:
                        gpos1 = eval(exonStarts[idx - 1]) - 1
                        gpos2 = eval(exonStarts[idx - 1])
                        gpos3 = eval(exonEnds[idx]) + 1
                        break
                    elif (cdsStart - eval(exonEnds[idx])) % 3 == 2:
                        gpos1 = eval(exonStarts[idx - 1])
                        gpos2 = eval(exonEnds[idx]) + 1
                        gpos3 = eval(exonEnds[idx]) + 2
                        break
                else:
                    gpos1 = cdsStart - pcPos * 3 + 1
                    gpos2 = cdsStart - pcPos * 3 + 2
                    gpos3 = cdsStart - pcPos * 3 + 3
                    break
            elif idx > edx:
                if pcPos - (eval(exonStarts[idx]) - eval(exonEnds[idx]) + eval(exonFrames[idx])) // 3 > 1:
                    pcPos -= (eval(exonStarts[idx]) - eval(exonEnds[idx]) + eval(exonFrames[idx])) // 3
                elif pcPos - (
                        eval(exonStarts[idx]) - eval(exonEnds[idx]) + eval(exonFrames[idx])) // 3 == 1 and (
                        eval(exonStarts[idx]) - eval(exonEnds[idx]) + eval(exonFrames[idx])) % 3 != 0:
                    if (eval(exonStarts[idx]) - eval(exonEnds[idx]) + eval(exonFrames[idx])) % 3 == 1:
                        gpos1 = eval(exonStarts[idx - 1]) - 1
                        gpos2 = eval(exonStarts[idx - 1])
                        gpos3 = eval(exonEnds[idx]) + 1
                        break
                    elif (eval(exonStarts[idx]) - eval(exonEnds[idx]) + eval(exonFrames[idx])) % 3 == 2:
                        gpos1 = eval(exonStarts[idx - 1])
                        gpos2 = eval(exonEnds[idx]) + 1
                        gpos3 = eval(exonEnds[idx]) + 2
                        break
                else:
                    gpos1 = eval(exonStarts[idx]) + eval(exonFrames[idx]) - pcPos * 3 + 1
                    gpos2 = eval(exonStarts[idx]) + eval(exonFrames[idx]) - pcPos * 3 + 2
                    gpos3 = eval(exonStarts[idx]) + eval(exonFrames[idx]) - pcPos * 3 + 3
                    break
            elif idx == edx:
                if pcPos - (eval(exonStarts[idx]) - cdsEns + (eval(exonFrames[idx])) % 3) // 3 > 1 or \
                        (pcPos - (eval(exonStarts[idx]) - cdsEns + eval(exonFrames[idx])) // 3 == 1 and (
                                eval(exonStarts[idx]) - cdsEns + eval(exonFrames[idx])) % 3 != 0):
                    print(gene, '超出转录本长度')
                    gpos1, gpos2, gpos3 = 0, 0, 0
                else:
                    gpos1 = eval(exonStarts[idx]) + eval(exonFrames[idx]) - pcPos * 3 + 1
                    gpos2 = eval(exonStarts[idx]) + eval(exonFrames[idx]) - pcPos * 3 + 2
                    gpos3 = eval(exonStarts[idx]) + eval(exonFrames[idx]) - pcPos * 3 + 3
                    break

    # print(gpos1, gpos2, gpos3, alt1, alt2, alt3)
    varlist = []
    if gpos1 and (alt1 or alt2 or alt3):
        # print(chrom + '\t' + str(gpos1))
        # print(chrom + '\t' + str(gpos2))
        # print(chrom + '\t' + str(gpos3))
        for i in range(len(alt1)):
            varlist.append(chrom + '\t' + str(gpos1) + '\t' + str(gpos1) + '\t' + ref1 + '\t' + alt1[i])
        for i in range(len(alt2)):
            varlist.append(chrom + '\t' + str(gpos2) + '\t' + str(gpos2) + '\t' + ref2 + '\t' + alt2[i])
        for i in range(len(alt3)):
            varlist.append(chrom + '\t' + str(gpos3) + '\t' + str(gpos3) + '\t' + ref3 + '\t' + alt3[i])
    return varlist
