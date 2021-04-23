na = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

aa1to3 = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
          'G': 'GLY', 'H': 'HIS', 'K': 'LYS', 'I': 'ILE', 'L': 'LEU',
          'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
          'S': 'SER', 'T': 'THR', 'V': 'VAL', 'Y': 'TYR', 'W': 'TRP', '*': 'Ter'}

transl_table_1 = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                  'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                  'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
                  'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
                  'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                  'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                  'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                  'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                  'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                  'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                  'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
                  'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
                  'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                  'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                  'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                  'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}
transl_table_5 = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                  'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                  'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': 'W',
                  'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
                  'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                  'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                  'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                  'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                  'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                  'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                  'ATA': 'M', 'ACA': 'T', 'AAA': 'K', 'AGA': 'S',
                  'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'S',
                  'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                  'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                  'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                  'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}


def format_input(dna):
    dna = dna.upper().replace('U', 'T')
    return dna


def tran(dna, transl_table):
    tran_aa = ['', '', '']
    for j in range(3):
        for i in range(int(len(dna) / 3)):
            if len(dna) - 3 * i > 2 + j:
                tran_aa[j] += (transl_table[dna[i * 3 + j:i * 3 + 3 + j]])
    return tran_aa


def score_tran_aa(tran_aa):
    score = []
    if tran_aa[0] == 'M':
        if_orf = [1]
        star_orf = 1
        end_orf = 0
    else:
        if_orf = [0]
        star_orf = 0
        end_orf = 1
    count = 0
    for aa in tran_aa:
        count += 1
        if aa == 'M' and star_orf == 0 and end_orf == 1:
            if_orf.append(1)
            score.append(count)
            count = 0
            end_orf = 0
            star_orf = 1
        if aa == '*' and star_orf == 1 and end_orf == 0:
            if_orf.append(0)
            score.append(count)
            count = 0
            end_orf = 1
            star_orf = 0
    score.append(count)
    print(score)
    print(if_orf)
    return score


def main():
    dna = 'atgataaaaataactttatctatattaattttttttatatcaattataatatatttt' \
          'ataaatcacccactatcaataggaatattaattttacttcaaacattattaacatgttta' \
          'ttgtcagggatattaatcaaaacttattgattttcatatattcttttcctagttttttta' \
          'ggaggattattagttttatttatttatgtatcaagtattgcatcaaatgaaatattttta' \
          'ttatcaaataatataaaaatcacattcattgtaatattaattataataatcagatttcaa' \
          'tttatatttactaaaaatttaaattgaataaatttaattaataactcagaaataaataat' \
          'tttttaaattttatatttttcaataatgaaaataaaattaatttaaataaactatataat' \
          'aataactcatctatattaatattattattaattatttatttatttatcacattaattgca' \
          'gttgtaaaaatcactaatattttttttggacctttacgaacttttaataattaa'
    dna = format_input(dna)
    dna_t = dna[::-1].translate(str.maketrans(na))

    tran_aa_5to3 = tran(dna, transl_table_5)
    tran_aa_3to5 = tran(dna_t, transl_table_5)
    for aa in tran_aa_5to3 + tran_aa_3to5:
        print(aa)

    test = tran_aa_5to3 + tran_aa_3to5
    score = score_tran_aa(test[0])


if __name__ == '__main__':
    main()
