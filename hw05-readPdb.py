aa3_1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
         'GLY': 'G', 'HIS': 'H', 'LYS': 'K', 'ILE': 'I', 'LEU': 'L',
         'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
         'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TYR': 'Y', 'TRP': 'W'}

filename = '6xre'
seq = ''
lastL = 'A'
lastN = 0
with open(filename + '.pdb', 'r', encoding='utf-8') as fr:
    with open(filename + '.fasta', 'w', encoding='utf-8') as fw:
        fw.write('> ' + filename + ' | Chain A | ' + '\n')
        for line in fr:
            if line[0:4] == 'ATOM':
                col = line.split()
                if not col[5] == lastN:
                    seq = seq + aa3_1[col[3]]
                    if not col[4][0] == lastL:
                        fw.write(seq + '\n' + '> ' + filename + ' | Chain ' + col[4][0] + ' |\n')
                        seq = ''
                lastL = col[4][0]
                lastN = col[5]
        fw.write(seq + '\n')
