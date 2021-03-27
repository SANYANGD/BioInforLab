
# Needleman-Wunsch 算法
# 区分开空位，延伸空位

import numpy as np

a_n = {'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4, 'Q': 5, 'E': 6, 'G': 7, 'H': 8, 'I': 9, 'L': 10,
       'K': 11, 'M': 12, 'F': 13, 'P': 14, 'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19, 'B': 20,
       'Z': 21, 'X': 22, '*': 23}
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
gap_penalty_e = -1


def nw(str1, str2):
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


def nwPlus(str1, str2):
    len_str1 = len(str1) + 1
    len_str2 = len(str2) + 1
    newstr1 = newstr2 = ''
    matrix_score = np.zeros([len_str2, len_str1])
    matrix_path = np.zeros([len_str2, len_str1])
    for m in range(len_str1):
        if m is 1:
            matrix_score[0, m] = gap_penalty
        elif m is not 1 and m is not 0:
            matrix_score[0, m] = matrix_score[0, m-1] + gap_penalty_e
    for n in range(len_str2):
        if n is 1:
            matrix_score[n, 0] = gap_penalty
        elif n is not 1 and n is not 0:
            matrix_score[n, 0] = matrix_score[n-1, 0] + gap_penalty_e

    for i in range(1, len_str1):
        for j in range(1, len_str2):
            t1 = a_n[str1[i - 1]]
            t2 = a_n[str2[j - 1]]

            if matrix_path[j, i-1] == 4:
                a = matrix_score[j, i - 1] + gap_penalty_e
            else:
                a = matrix_score[j, i - 1] + gap_penalty
            if matrix_path[j - 1, i] == 2:
                b = matrix_score[j - 1, i] + gap_penalty_e
            else:
                b = matrix_score[j - 1, i] + gap_penalty
            c = matrix_score[j - 1, i - 1] + blosum[t1, t2]
            matrix_score[j, i] = max(a, b, c)

            # from left up xie
            # 100 4; 010 2; 001 1; 101 5; 011 3
            path = 0b000
            if max(a, b, c) is a:
                path = path | 0b100
            elif max(a, b, c) is b:
                path = path | 0b010
            # elif max(a, b, c) is c:
            else:
                path = path | 0b001
            matrix_path[j, i] = path

    print(matrix_score)
    print(matrix_path)

    m = len_str2 - 1
    n = len_str1 - 1
    while m != 0 or n != 0:
        if matrix_path[m, n] == 4:
            newstr1 = str1[n - 1].__add__(newstr1)
            newstr2 = '*'.__add__(newstr2)
            n = n - 1
        elif matrix_path[m, n] == 2:
            newstr1 = '*'.__add__(newstr1)
            newstr2 = str2[m - 1].__add__(newstr2)
            m = m - 1
        elif matrix_path[m, n] == 1:
            newstr1 = str1[n - 1].__add__(newstr1)
            newstr2 = str2[m - 1].__add__(newstr2)
            m = m - 1
            n = n - 1
    print(newstr1)
    print(newstr2)


def main():
    #str1 = 'IPGAWD'
    #str2 = 'VGAWAD'

    str1 = 'MSRQFTCKSGAAAKGGFSGCSAVLSGGSSSSFRAGSKGLSGGFGSRSLYSLGGVRSLNVASGSGKSGGYGFGRGRASGFAGSMFGSVALGPVCPTVCPPGGIHQVTVNESLLAPLNVELDPEIQKVRAQEREQIKALNNKFASFIDKVRFLEQQNQVLETKWELLQQLDLNNCKNNLEPILEGYISNLRKQLETLSGDRVRLDSELRNVRDVVEDYKKRYEEEINKRTAAENEFVLLKKDVDAAYANKVELQAKVESMDQEIKFFRCLFEAEITQIQSHISDMSVILSMDNNRNLDLDSIIDEVRTQYEEIALKSKAEAEALYQTKFQELQLAAGRHGDDLKNTKNEISELTRLIQRIRSEIENVKKQASNLETAIADAEQRGDNALKDARAKLDELEGALHQAKEELARMLREYQELMSLKLALDMEIATYRKLLESEECRMSGEFPSPVSISIISSTSGGSVYGFRPSMVSGGYVANSSNCISGVCSVRGGEGRSRGSANDYKDTLGKGSSLSAPSKKTSR'
    str2 = 'MSRQFTCKSGASNRGFSGCSAVLSGGSSSSYRAGGKGLSGGFGSRSLYSLGGGRSITLNMASGSGKNGGFGFGRNRASGFAGSIFGSVALGPVCPAVCPPGGIHQVTVNESLLAPLNVELDPEIQKVRAQEREQIKALNNKFASFIDKVRFLEQQNQVLQTKWELLQQLDLNNCKNNLEPILEGHISNMRKQLETLSGDRVRLDSELRNVRDVVEDYKKKYEEEINRRTAAENEFVLLKKDVDAAYANKVELQAKVDTMDQDIKFFKCLFEAEMAQIQSHISDMSVILSMDNNRNLDLDSIIDEVRAQYEEIALKSKAEAEALYQTKFQELQLAAGRHGDDLKNTKNEITELTRFIQRLRSEIENAKKQASNLETAIADAEQRGDSALKDARAKLDELEGALHQAKEELARMLREYQELMSLKLALDMEIATYRKLLESEECRMSGEYSSPVSISIISSTSGSGGYGFRPSTVSGGYVANSTSCISGVCSVRGGENRSRGSASDYKDTLTKGSSLSTPSKKGGR'

    # nw(str1, str2)
    nwPlus(str1, str2)


if __name__ == '__main__':
    main()
