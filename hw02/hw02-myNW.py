
# Needleman-Wunsch 算法

import numpy as np

a_n = {'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4, 'Q': 5, 'E': 6, 'G': 7, 'H': 8, 'I': 9, 'L': 10,
       'K': 11, 'M': 12, 'F': 13, 'P': 14, 'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19, 'B': 20,
       'Z': 21, 'X': 22, '*': 23}
# #  Matrix made by matblas from blosum62.iij
# #  * column uses minimum score
# #  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
# #  Blocks Database = /data/blocks_5.0/blocks.dat
# #  Cluster Percentage: >= 62
# #  Entropy =   0.6979, Expected =  -0.5209
#    A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
# A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
# R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
# N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
# D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
# C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
# Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
# E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
# G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
# H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
# I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
# L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
# K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
# M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
# F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
# P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
# S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
# T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
# W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
# Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
# V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
# B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
# Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
# X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
# * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
blosum_62 = [[4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4],
             [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4],
             [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4],
             [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4],
             [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4],
             [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4],
             [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4],
             [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4],
             [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4],
             [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4],
             [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4],
             [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4],
             [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4],
             [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4],
             [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4],
             [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4],
             [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4],
             [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4],
             [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4],
             [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4],
             [-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4],
             [-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4],
             [0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4],
             [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1]]
blosum_45 = [[5, -2, -1, -2, -1, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -2, -2, 0, -1, -1, -1, -1, -5],
             [-2, 7, 0, -1, -3, 1, 0, -2, 0, -3, -2, 3, -1, -2, -2, -1, -1, -2, -1, -2, -1, -3, 1, -1, -5],
             [-1, 0, 6, 2, -2, 0, 0, 0, 1, -2, -3, 0, -2, -2, -2, 1, 0, -4, -2, -3, 5, -3, 0, -1, -5],
             [-2, -1, 2, 7, -3, 0, 2, -1, 0, -4, -3, 0, -3, -4, -1, 0, -1, -4, -2, -3, 6, -3, 1, -1, -5],
             [-1, -3, -2, -3, 12, -3, -3, -3, -3, -3, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -2, -2, -3, -1, -5],
             [-1, 1, 0, 0, -3, 6, 2, -2, 1, -2, -2, 1, 0, -4, -1, 0, -1, -2, -1, -3, 0, -2, 4, -1, -5],
             [-1, 0, 0, 2, -3, 2, 6, -2, 0, -3, -2, 1, -2, -3, 0, 0, -1, -3, -2, -3, 1, -3, 5, -1, -5],
             [0, -2, 0, -1, -3, -2, -2, 7, -2, -4, -3, -2, -2, -3, -2, 0, -2, -2, -3, -3, -1, -4, -2, -1, -5],
             [-2, 0, 1, 0, -3, 1, 0, -2, 10, -3, -2, -1, 0, -2, -2, -1, -2, -3, 2, -3, 0, -2, 0, -1, -5],
             [-1, -3, -2, -4, -3, -2, -3, -4, -3, 5, 2, -3, 2, 0, -2, -2, -1, -2, 0, 3, -3, 4, -3, -1, -5],
             [-1, -2, -3, -3, -2, -2, -2, -3, -2, 2, 5, -3, 2, 1, -3, -3, -1, -2, 0, 1, -3, 4, -2, -1, -5],
             [-1, 3, 0, 0, -3, 1, 1, -2, -1, -3, -3, 5, -1, -3, -1, -1, -1, -2, -1, -2, 0, -3, 1, -1, -5],
             [-1, -1, -2, -3, -2, 0, -2, -2, 0, 2, 2, -1, 6, 0, -2, -2, -1, -2, 0, 1, -2, 2, -1, -1, -5],
             [-2, -2, -2, -4, -2, -4, -3, -3, -2, 0, 1, -3, 0, 8, -3, -2, -1, 1, 3, 0, -3, 1, -3, -1, -5],
             [-1, -2, -2, -1, -4, -1, 0, -2, -2, -2, -3, -1, -2, -3, 9, -1, -1, -3, -3, -3, -2, -3, -1, -1, -5],
             [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -3, -1, -2, -2, -1, 4, 2, -4, -2, -1, 0, -2, 0, -1, -5],
             [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1, 2, 5, -3, -1, 0, 0, -1, -1, -1, -5],
             [-2, -2, -4, -4, -5, -2, -3, -2, -3, -2, -2, -2, -2, 1, -3, -4, -3, 15, 3, -3, -4, -2, -2, -1, -5],
             [-2, -1, -2, -2, -3, -1, -2, -3, 2, 0, 0, -1, 0, 3, -3, -2, -1, 3, 8, -1, -2, 0, -2, -1, -5],
             [0, -2, -3, -3, -1, -3, -3, -3, -3, 3, 1, -2, 1, 0, -3, -1, 0, -3, -1, 5, -3, 2, -3, -1, -5],
             [-1, -1, 5, 6, -2, 0, 1, -1, 0, -3, -3, 0, -2, -3, -2, 0, 0, -4, -2, -3, 5, -3, 1, -1, -5],
             [-1, -3, -3, -3, -2, -2, -3, -4, -2, 4, 4, -3, 2, 1, -3, -2, -1, -2, 0, 2, -3, 4, -2, -1, -5],
             [-1, 1, 0, 1, -3, 4, 5, -2, 0, -3, -2, 1, -1, -3, -1, 0, -1, -2, -2, -3, 1, -2, 5, -1, -5],
             [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5],
             [-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 1]]
blosum = np.array(blosum_62)
gap_penalty = -4


def dynamic_edit(str1, str2):
    len_str1 = len(str1) + 1
    len_str2 = len(str2) + 1
    newstr1 = newstr2 = ''
    matrix = np.zeros([len_str2, len_str1])
    matrix[0] = np.arange(0, gap_penalty * len_str1, gap_penalty)
    matrix[:, 0] = np.arange(0, gap_penalty * len_str2, gap_penalty)

    for i in range(1, len_str1):
        for j in range(1, len_str2):
            t1 = a_n[str1[i - 1]]
            t2 = a_n[str2[j - 1]]
            a = matrix[j, i - 1] + gap_penalty
            b = matrix[j - 1, i] + gap_penalty
            c = matrix[j - 1, i - 1] + blosum[t1, t2]
            matrix[j, i] = max(a, b, c)
    print(matrix)

    m = len_str2 - 1
    n = len_str1 - 1
    while m != 0 or n != 0:
        t1 = a_n[str1[n - 1]]
        t2 = a_n[str2[m - 1]]
        if matrix[m, n] == matrix[m, n - 1] + gap_penalty:
            newstr1 = str1[n - 1].__add__(newstr1)
            newstr2 = '*'.__add__(newstr2)
            n = n - 1
        elif matrix[m, n] == matrix[m - 1, n] + gap_penalty:
            newstr1 = '*'.__add__(newstr1)
            newstr2 = str2[m - 1].__add__(newstr2)
            m = m - 1
        elif matrix[m, n] == matrix[m - 1, n - 1] + blosum[t1, t2]:
            newstr1 = str1[n - 1].__add__(newstr1)
            newstr2 = str2[m - 1].__add__(newstr2)
            m = m - 1
            n = n - 1
    print(newstr1)
    print(newstr2)


def main():
    # str1 = 'IPGAWD'
    # str2 = 'VGAWAD'

    str1 = 'MSCQISCKSRGRGGGGGGFRGFSSGSAVVSGGSRRSTSSFSCLSRHGGGGGGFGGGGFGSRSLVGLGGTKSISISVAGGGGGFGAAGGFGGRGGGFGGGSSFGGGSGFSGGGFGGGGFGGGRFGGFGGPGGVGGLGGPGGFGPGGYPGGIHEVSVNQSLLQPLNVKVDPEIQNVKAQEREQIKTLNNKFASFIDKVRFLEQQNQVLQTKWELLQQMNVGTRPINLEPIFQGYIDSLKRYLDGLTAERTSQNSELNNMQDLVEDYKKKYEDEINKRTAAENDFVTLKKDVDNAYMIKVELQSKVDLLNQEIEFLKVLYDAEISQIHQSVTDTNVILSMDNSRNLDLDSIIAEVKAQYEEIAQRSKEEAEALYHSKYEELQVTVGRHGDSLKEIKIEISELNRVIQRLQGEIAHVKKQCKNVQDAIADAEQRGEHALKDARNKLNDLEEALQQAKEDLARLLRDYQELMNVKLALDVEIATYRKLLEGEECRMSGDLSSNVTVSVTSSTISSNVASKAAFGGSGGRGSSSGGGYSSGSSSYGSGGRQSGSRGGSGGGGSISGGGYGSGGGSGGRYGSGGGSKGGSISGGGYGSGGGKHSSGGGSRGGSSSGGGYGSGGGGSSSVKGSSGEAFGSSVTFSFR'
    str2 = 'MSCQISCRSRRGGGGGGGGGFRGFSSGSAVVSGGSRRSNTSFSCISRHGGGRGGSGGGGFGSQSLVGLGGYKSISSSVAGNSGGYGGSSFGGSSGFGGGRGFGGGQGFGGSGGFGGGSGFGGGQGFGGGSRFGGGSGFGGGGFGGGSFGGGRFGGGPGGFGGPGGFPGGGIHEVSVNQSLLQPLDVKVDPEIQNVKSQEREQIKTLNNKFASFIDKVRFLEQQNQVLRTKWELLQQLDVGSRTTNLDPIFQAYIGMLKKQVDRLSAERTSQESELNNMQDLVEDFKKKYEDEINKRTSAENDFVTIKKDVDSCYMDKTELQARLDILAQEVNFLRTLYDAELSQLQQDVTDTNVILSMDNNRNLDLDSIIAEVQNQYEMIAHKSKAESEELYHSKYEELQVTAVKHGDSLKEIKMEISELNRTIQRLQGEISHVKKQCKGVQDSIADAEQRGEHAIKDARGKLTDLEEALQQCREDLARLLRDYQELMNTKLSLDVEIATYRKLLEGEECRMSGDFSDNVSVSITSSTISSSVASKTGFGSGGQSSGGRGSYGGRGGGGGGGSTYGSGGRSSGSRGSGSGSGGGGYSSGGGSRGGSGGGYGSGGGSRGGSGGGYGSGGGSGSGGGYSSGGGSRGGSGGGGVSSGGGSRGGSSSGGGSRGGSSSGGGGYSSGGGSRGGSSSGGAGSSSEKGGSGSGEGCGSGVTFSFR'

    dynamic_edit(str1, str2)


if __name__ == '__main__':
    main()
