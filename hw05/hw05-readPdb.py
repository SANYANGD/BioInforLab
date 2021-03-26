aa3_1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
         'GLY': 'G', 'HIS': 'H', 'LYS': 'K', 'ILE': 'I', 'LEU': 'L',
         'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
         'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TYR': 'Y', 'TRP': 'W'}

# filename = '4a61'
filename = input('请输入PDB序号：')
seq = ''
lastL = 'A'
lastN = 0
changeM = 0
with open(filename + '.pdb', 'r', encoding='utf-8') as fr:
    with open(filename + '.fasta', 'w', encoding='utf-8') as fw:
        fw.write('> ' + filename + ' | Chain A | ' + '\n')
        for line in fr:
            if line[0:4] == 'ATOM':
                col = line.split()
                if not col[5] == lastN:  # 检测编号是否改变
                    seq = seq + aa3_1[col[3]]

                    if not col[4][0] == lastL:  # 检测亚基是否改变，亚基变化，开新seq
                        fw.write(seq + '\n' + '> ' + filename + ' | Chain ' + col[4][0] + ' |\n')
                        seq = ''
                        changeM = 0  # 编号从新计数了

                lastL = col[4][0]  # 当编号大于999时，会出现如 A1002 这种问题

                if changeM == 0:
                    lastN = col[5]
                    if int(col[5]) >= 999:
                        changeM = 1  # 此时编号出现如 A1002 这种问题，换方法提区编号
                elif changeM == 1:
                    lastN = col[4][1:]

        fw.write(seq + '\n')
