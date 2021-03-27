
amino = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
         'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']
blosum = [[4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4],
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
gap_penalty = -11
gap_penalty_e = -1


def find_score(str1, str2, i, j):
    index1 = amino.index(str2[i - 1])
    index2 = amino.index(str1[j - 1])
    return blosum[index1][index2]


def creatx(i, j):
    if i > 0 and j == 0:
        return -float("inf")
    else:
        if j > 0:
            return gap_penalty + gap_penalty_e * j
        else:
            return 0


def creaty(i, j):
    if j > 0 and i == 0:
        return -float("inf")
    else:
        if i > 0:
            return gap_penalty + gap_penalty_e * i
        else:
            return 0


def creatm(i, j):
    if j == 0 and i == 0:
        return 0
    else:
        if j == 0 or i == 0:
            return -float("inf")
        else:
            return 0


# 生成距离矩阵
def score_matrix(str1, str2):
    len_str2 = len(str2) + 1
    len_str1 = len(str1) + 1
    X = [[creatx(i, j) for j in range(0, len_str1)] for i in range(0, len_str2)]
    Y = [[creaty(i, j) for j in range(0, len_str1)] for i in range(0, len_str2)]
    M = [[creatm(i, j) for j in range(0, len_str1)] for i in range(0, len_str2)]
    for j in range(1, len_str1):
        for i in range(1, len_str2):
            # 因为同状态之间不能相互转换（X不能转向Y)，所以X有两种情况存在。X代表竖直情况
            # 第一种是X转向X，即X[i][j]=X[i][j-1]+e
            # 第二种是M态转向，X[x][j-1]=M[i][j-1]+gap_penalty+gap_penalty_e
            X[i][j] = max((gap_penalty + gap_penalty_e + M[i][j - 1]), (gap_penalty_e + X[i][j - 1]))
            # 和上述X状态一致，实则X、Y为相同状态。只不过Y[i][j]中是另一个平面的情况。Y代表水平情况
            # 两者只是简单i,j的替换
            Y[i][j] = max((gap_penalty + gap_penalty_e + M[i - 1][j]), (gap_penalty_e + Y[i - 1][j]))
            # M状态则存在三种情况
            # 一个是自身的转换，即M态到M态，这个时候就是对角线的转换。加上上述中的mismatch或则match积分
            # 第二种和第三种则是对于X、Y状态的转换。这种转换只能发生在当X、Y为边缘时，即为M本身的情况
            M[i][j] = max(find_score(str1, str2, i, j) + M[i - 1][j - 1], X[i][j], Y[i][j])
    return X, Y, M


def back(str1, str2, X, Y, M):
    seq1 = ''
    seq2 = ''
    len_str2 = len(str2)
    len_str1 = len(str1)
    while len_str2 > 0 or len_str1 > 0:
        # 当i，j不为0，即不在边界上时。且M状态是直接转换M态。说明对角线进行转换。蛋白质中的氨基酸序列中相匹配（不一定是相等）,但是为最优。
        if len_str2 > 0 and len_str1 > 0 and \
                M[len_str2][len_str1] == M[len_str2 - 1][len_str1 - 1] + find_score(str1, str2, len_str2, len_str1):
            seq1 += str1[len_str1 - 1]
            seq2 += str2[len_str2 - 1]
            len_str2 -= 1
            len_str1 -= 1
        # 在竖直状态下为空位
        elif len_str2 > 0 and M[len_str2][len_str1] == Y[len_str2][len_str1]:
            # 当i>0,且j=0时，我们只对i进行相减.打分矩阵与Y矩阵值一致，说明为存在空位
            seq1 += '-'
            seq2 += str2[len_str2 - 1]
            len_str2 -= 1
        # 在水平状态下为空位
        elif len_str1 > 0 and M[len_str2][len_str1] == X[len_str2][len_str1]:
            seq1 += str1[len_str1 - 1]
            seq2 += '-'
            len_str1 -= 1

    seq1_f = ''.join([seq1[j] for j in range(-1, -(len(seq1) + 1), -1)])
    seq2_f = ''.join([seq2[j] for j in range(-1, -(len(seq2) + 1), -1)])

    return seq1_f, seq2_f


def main():
    #str1 = 'IPGAWD'
    #str2 = 'VGAWAD'
    str1 = 'MSCQISCKSRGRGGGGGGFRGFSSGSAVVSGGSRRSTSSFSCLSRHGGGGGGFGGGGFGSRSLVGLGGTKSISISVAGGGGGFGAAGGFGGRGGGFGGGSSFGGGSGFSGGGFGGGGFGGGRFGGFGGPGGVGGLGGPGGFGPGGYPGGIHEVSVNQSLLQPLNVKVDPEIQNVKAQEREQIKTLNNKFASFIDKVRFLEQQNQVLQTKWELLQQMNVGTRPINLEPIFQGYIDSLKRYLDGLTAERTSQNSELNNMQDLVEDYKKKYEDEINKRTAAENDFVTLKKDVDNAYMIKVELQSKVDLLNQEIEFLKVLYDAEISQIHQSVTDTNVILSMDNSRNLDLDSIIAEVKAQYEEIAQRSKEEAEALYHSKYEELQVTVGRHGDSLKEIKIEISELNRVIQRLQGEIAHVKKQCKNVQDAIADAEQRGEHALKDARNKLNDLEEALQQAKEDLARLLRDYQELMNVKLALDVEIATYRKLLEGEECRMSGDLSSNVTVSVTSSTISSNVASKAAFGGSGGRGSSSGGGYSSGSSSYGSGGRQSGSRGGSGGGGSISGGGYGSGGGSGGRYGSGGGSKGGSISGGGYGSGGGKHSSGGGSRGGSSSGGGYGSGGGGSSSVKGSSGEAFGSSVTFSFR'
    str2 = 'MSCQISCRSRRGGGGGGGGGFRGFSSGSAVVSGGSRRSNTSFSCISRHGGGRGGSGGGGFGSQSLVGLGGYKSISSSVAGNSGGYGGSSFGGSSGFGGGRGFGGGQGFGGSGGFGGGSGFGGGQGFGGGSRFGGGSGFGGGGFGGGSFGGGRFGGGPGGFGGPGGFPGGGIHEVSVNQSLLQPLDVKVDPEIQNVKSQEREQIKTLNNKFASFIDKVRFLEQQNQVLRTKWELLQQLDVGSRTTNLDPIFQAYIGMLKKQVDRLSAERTSQESELNNMQDLVEDFKKKYEDEINKRTSAENDFVTIKKDVDSCYMDKTELQARLDILAQEVNFLRTLYDAELSQLQQDVTDTNVILSMDNNRNLDLDSIIAEVQNQYEMIAHKSKAESEELYHSKYEELQVTAVKHGDSLKEIKMEISELNRTIQRLQGEISHVKKQCKGVQDSIADAEQRGEHAIKDARGKLTDLEEALQQCREDLARLLRDYQELMNTKLSLDVEIATYRKLLEGEECRMSGDFSDNVSVSITSSTISSSVASKTGFGSGGQSSGGRGSYGGRGGGGGGGSTYGSGGRSSGSRGSGSGSGGGGYSSGGGSRGGSGGGYGSGGGSRGGSGGGYGSGGGSGSGGGYSSGGGSRGGSGGGGVSSGGGSRGGSSSGGGSRGGSSSGGGGYSSGGGSRGGSSSGGAGSSSEKGGSGSGEGCGSGVTFSFR'
    X, Y, M = score_matrix(str1, str2)
    seq1, seq2 = back(str1, str2, X, Y, M)
    print("Alignment Score:" + str(M[len(str2)][len(str1)]))
    print("seq1:" + seq1)
    print("seq2:" + seq2)


if __name__ == '__main__':
    main()
