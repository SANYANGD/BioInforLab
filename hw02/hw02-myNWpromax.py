gap_penalty = -11  # gap opening penalty
gap_penalty_e = -1  # gap extending penalty
MIN = -float("inf")  # 用来进行边值初始化
amino = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z',
         'X', '*']
blosum = [
    [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4],
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
    [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1],
]


# return match or mismatch score
def _match(s, t, i, j):  # 该函数用于对于双序列S、T中不同氨基酸的match或则mismatch打分
    # 读取t字符串中t[i-1]字符对于amino对于的顺序值
    index1 = amino.index(t[i - 1])
    # 读取s字符串中t[j-1]字符对于amino对于的顺序值
    index2 = amino.index(s[j - 1])
    # index2=____________________
    return blosum[index1][index2]


# initializers for matrices
# 对于X状态情况下初始化，gap  extending=d+e*(n-1)  同等状态下的情况
def _init_x(i, j):
    if i > 0 and j == 0:
        return MIN
    else:
        if j > 0:
            return gap_penalty + gap_penalty_e * j
        else:
            return 0


# 对于Y状态情况下初始化，gap  extending=d+e*(n-1)  同等状态下的情况
def _init_y(i, j):
    if j > 0 and i == 0:
        return MIN
    else:
        if i > 0:
            return gap_penalty + gap_penalty_e * i
        else:
            return 0


# 对M状态进行初始化
def _init_m(i, j):
    if j == 0 and i == 0:
        return 0
    else:
        if j == 0 or i == 0:
            return MIN
        else:
            return 0


# 生成距离矩阵
def distance_matrix(s, t):
    dim_i = len(t) + 1
    dim_j = len(s) + 1
    # abuse list comprehensions to create matrices
    # 通过创建X、Y、M矩阵来实现创建三维空间
    X = [[_init_x(i, j) for j in range(0, dim_j)] for i in range(0, dim_i)]
    # print(X)
    Y = [[_init_y(i, j) for j in range(0, dim_j)] for i in range(0, dim_i)]
    M = [[_init_m(i, j) for j in range(0, dim_j)] for i in range(0, dim_i)]
    # print(M)
    for j in range(1, dim_j):
        for i in range(1, dim_i):
            # 因为同状态之间不能相互转换（X不能转向Y)，所以X有两种情况存在。X代表竖直情况
            # 第一种是X转向X，即X[i][j]=X[i][j-1]+e
            # 第二种是M态转向，X[x][j-1]=M[i][j-1]+gap_penalty+gap_penalty_e
            X[i][j] = max((gap_penalty + gap_penalty_e + M[i][j - 1]), (gap_penalty_e + X[i][j - 1]))
            # 和上述X状态一致，实则X、Y为相同状态。只不过Y[i][j]中是另一个平面的情况。Y代表水平情况
            # 两者只是简单i,j的替换
            # Y[i][j] = max(___________________, (gap_penalty_e + Y[i-1][j]))
            Y[i][j] = max((gap_penalty + gap_penalty_e + M[i - 1][j]), (gap_penalty_e + Y[i - 1][j]))
            # M状态则存在三种情况
            # 一个是自身的转换，即M态到M态，这个时候就是对角线的转换。加上上述中的mismatch或则match积分
            # 第二种和第三种则是对于X、Y状态的转换。这种转换只能发生在当X、Y为边缘时，即为M本身的情况
            # M[i][j] = max(_match(s, t, i, j) + M[i-1][j-1], ___________, ____________)
            # M[i][j] = max(_match(s, t, i, j) + M[i - 1][j - 1],X[i-1][j-1]+_match(s, t, i, j),Y[i-1][j-1]+_match(s, t, i, j))
            M[i][j] = max(_match(s, t, i, j) + M[i - 1][j - 1], X[i][j],
                          Y[i][j])
    return [X, Y, M]


# 回溯，得到双序列匹配的字符串
def backtrace(s, t, X, Y, M):
    sequ1 = ''
    sequ2 = ''
    i = len(t)
    j = len(s)
    while (i > 0 or j > 0):
        # 当i，j不为0，即不在边界上时。且M状态是直接转换M态。说明对角线进行转换。
        # 蛋白质中的氨基酸序列中相匹配（不一定是相等）,但是为最优。
        if (i > 0 and j > 0 and M[i][j] == M[i - 1][j - 1] + _match(s, t, i, j)):
            sequ1 += s[j - 1]
            sequ2 += t[i - 1]
            i -= 1
            j -= 1
        # 当j=0且i不为0时，说明只可以进行水平转移，说明在竖直状态下为空位
        elif (i > 0 and M[i][j] == Y[i][j]):
            # sequ1 += ______
            # 当i>0,且j=0时，我们只对i进行相减.打分矩阵与Y矩阵值一致，说明为存在空位
            sequ1 += '_'
            sequ2 += t[i - 1]
            i -= 1
        # 当i=0且j不为0时，说明只可以进行竖直转移，说明在水平状态下为空位
        elif (j > 0 and M[i][j] == X[i][j]):
            # sequ1 += _______
            # 当j>0,且i=0时，我们只对j进行相减，打分矩阵与X矩阵值一致，说明为存在空位
            sequ1 += s[j - 1]
            sequ2 += '_'
            j -= 1
    # 通过回溯后，倒序输出
    sequ1r = ''.join([sequ1[j] for j in range(-1, -(len(sequ1) + 1), -1)])
    sequ2r = ''.join([sequ2[j] for j in range(-1, -(len(sequ2) + 1), -1)])

    return [sequ1r, sequ2r]


def main():
    seq1 = 'MSRQFTCKSGAAAKGGFSGCSAVLSGGSSSSFRAGSKGLSGGFGSRSLYSLGGVRSLNVASGSGKSGGYGFGRGRASGFAGSMFGSVALGPVCPTVCPPGGIHQVTVNESLLAPLNVELDPEIQKVRAQEREQIKALNNKFASFIDKVRFLEQQNQVLETKWELLQQLDLNNCKNNLEPILEGYISNLRKQLETLSGDRVRLDSELRNVRDVVEDYKKRYEEEINKRTAAENEFVLLKKDVDAAYANKVELQAKVESMDQEIKFFRCLFEAEITQIQSHISDMSVILSMDNNRNLDLDSIIDEVRTQYEEIALKSKAEAEALYQTKFQELQLAAGRHGDDLKNTKNEISELTRLIQRIRSEIENVKKQASNLETAIADAEQRGDNALKDARAKLDELEGALHQAKEELARMLREYQELMSLKLALDMEIATYRKLLESEECRMSGEFPSPVSISIISSTSGGSVYGFRPSMVSGGYVANSSNCISGVCSVRGGEGRSRGSANDYKDTLGKGSSLSAPSKKTSR'
    seq2 = 'MSRQFTCKSGASNRGFSGCSAVLSGGSSSSYRAGGKGLSGGFGSRSLYSLGGGRSITLNMASGSGKNGGFGFGRNRASGFAGSIFGSVALGPVCPAVCPPGGIHQVTVNESLLAPLNVELDPEIQKVRAQEREQIKALNNKFASFIDKVRFLEQQNQVLQTKWELLQQLDLNNCKNNLEPILEGHISNMRKQLETLSGDRVRLDSELRNVRDVVEDYKKKYEEEINRRTAAENEFVLLKKDVDAAYANKVELQAKVDTMDQDIKFFKCLFEAEMAQIQSHISDMSVILSMDNNRNLDLDSIIDEVRAQYEEIALKSKAEAEALYQTKFQELQLAAGRHGDDLKNTKNEITELTRFIQRLRSEIENAKKQASNLETAIADAEQRGDSALKDARAKLDELEGALHQAKEELARMLREYQELMSLKLALDMEIATYRKLLESEECRMSGEYSSPVSISIISSTSGSGGYGFRPSTVSGGYVANSTSCISGVCSVRGGENRSRGSASDYKDTLTKGSSLSTPSKKGGR'

    # 需要将输入字符串转换为3个状态
    [X, Y, M] = distance_matrix(seq1, seq2)
    print(X)
    print(M)
    [str1, str2] = backtrace(seq1, seq2, X, Y, M)

    print("Alignment Score:" + str(M[len(seq2)][len(seq1)]))
    print(str1)
    print(str2)


if __name__ == '__main__':
    main()
