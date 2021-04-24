import re
import tkinter as tk
import numpy as np

# 输入窗体
window = tk.Tk()  # 建立窗口window
window.title('DNA/RAN Translate')  # 窗口名称

# 碱基配对表
na = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

# 氨基酸单字母-三字母对应表
aa1to3 = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
          'G': 'GLY', 'H': 'HIS', 'K': 'LYS', 'I': 'ILE', 'L': 'LEU',
          'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
          'S': 'SER', 'T': 'THR', 'V': 'VAL', 'Y': 'TYR', 'W': 'TRP', '*': 'Ter'}

# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG1
# 根据网站得到的多种transl_table
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
transl_table_2 = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                  'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                  'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': 'W',
                  'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
                  'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                  'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                  'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                  'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                  'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                  'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                  'ATA': 'M', 'ACA': 'T', 'AAA': 'K', 'AGA': '*',
                  'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': '*',
                  'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                  'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                  'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                  'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}
transl_table_3 = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                  'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                  'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': 'W',
                  'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
                  'CTT': 'T', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                  'CTC': 'T', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                  'CTA': 'T', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                  'CTG': 'T', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                  'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                  'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                  'ATA': 'M', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
                  'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
                  'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                  'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                  'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                  'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}
transl_table_4 = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
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
transl_table_6 = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                  'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                  'TTA': 'L', 'TCA': 'S', 'TAA': 'Q', 'TGA': '*',
                  'TTG': 'L', 'TCG': 'S', 'TAG': 'Q', 'TGG': 'W',
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
transl_table_9 = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                  'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                  'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': 'W',
                  'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
                  'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                  'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                  'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                  'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                  'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                  'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                  'ATA': 'I', 'ACA': 'T', 'AAA': 'N', 'AGA': 'S',
                  'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'S',
                  'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                  'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                  'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                  'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}
transl_table_10 = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                   'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                   'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': 'C',
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
transl_table_12 = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                   'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                   'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
                   'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
                   'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                   'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                   'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                   'CTG': 'S', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                   'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                   'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                   'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
                   'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
                   'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                   'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                   'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                   'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}
transl_table_13 = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                   'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                   'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': 'W',
                   'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
                   'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                   'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                   'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                   'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                   'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                   'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                   'ATA': 'M', 'ACA': 'T', 'AAA': 'K', 'AGA': 'G',
                   'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'G',
                   'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                   'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                   'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                   'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}
transl_table_14 = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                   'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                   'TTA': 'L', 'TCA': 'S', 'TAA': 'Y', 'TGA': 'W',
                   'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
                   'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                   'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                   'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                   'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                   'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                   'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                   'ATA': 'I', 'ACA': 'T', 'AAA': 'N', 'AGA': 'S',
                   'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'S',
                   'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                   'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                   'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                   'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}
transl_table_16 = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                   'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                   'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
                   'TTG': 'L', 'TCG': 'S', 'TAG': 'L', 'TGG': 'W',
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
transl_table_21 = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                   'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                   'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': 'W',
                   'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
                   'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                   'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                   'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                   'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                   'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                   'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                   'ATA': 'M', 'ACA': 'T', 'AAA': 'N', 'AGA': 'S',
                   'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'S',
                   'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                   'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                   'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                   'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}
transl_table_22 = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                   'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                   'TTA': 'L', 'TCA': '*', 'TAA': '*', 'TGA': '*',
                   'TTG': 'L', 'TCG': 'S', 'TAG': 'L', 'TGG': 'W',
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
transl_table = [transl_table_1, transl_table_2, transl_table_3, transl_table_4, transl_table_5,
                transl_table_6, transl_table_9, transl_table_10, transl_table_12, transl_table_13,
                transl_table_14, transl_table_16, transl_table_21, transl_table_22]


# 格式化DNA输入
# def format_input():
#     # dna = 'atgataaaaataactttatctatattaattttttttatatcaattataatatatttt' \
#     #       'ataaatcacccactatcaataggaatattaattttacttcaaacattattaacatgttta' \
#     #       'ttgtcagggatattaatcaaaacttattgattttcatatattcttttcctagttttttta' \
#     #       'ggaggattattagttttatttatttatgtatcaagtattgcatcaaatgaaatattttta' \
#     #       'ttatcaaataatataaaaatcacattcattgtaatattaattataataatcagatttcaa' \
#     #       'tttatatttactaaaaatttaaattgaataaatttaattaataactcagaaataaataat' \
#     #       'tttttaaattttatatttttcaataatgaaaataaaattaatttaaataaactatataat' \
#     #       'aataactcatctatattaatattattattaattatttatttatttatcacattaattgca' \
#     #       'gttgtaaaaatcactaatattttttttggacctttacgaacttttaataattaa'
#     dna = input('Please enter a protokaryotic DNA or RNA sequence - numbers and blanks are ignored:\n')
#     dna = dna.upper().replace('U', 'T')
#     dna = re.sub('[^ATGC]', '', dna)
#     return dna


# 结果输出
def format_output(result):
    print('\nResults of translation (Ranked by recommendation):')
    print('%-10s %-s' % ('score', 'aa'))
    for r in result:
        print('%-10s %-s' % (r[1], r[0]))


# DNA to AA
def tran(dna, transl_table):
    tran_aa = ['', '', '']
    for j in range(3):
        for i in range(int(len(dna) / 3)):
            if len(dna) - 3 * i > 2 + j:
                tran_aa[j] += (transl_table[dna[i * 3 + j:i * 3 + 3 + j]])
    return tran_aa


# 计算得到的AA序列结构
def num_tran_aa(tran_aa):
    num = []
    if tran_aa[0] == 'M':
        if_orf = [1]  # orf判定数组
        star_orf = 1
        end_orf = 0
    else:
        if_orf = [-1]
        star_orf = 0
        end_orf = 1
    count = 0
    for aa in tran_aa:
        count += 1
        if aa == 'M' and star_orf == 0 and end_orf == 1:
            if_orf.append(1)
            num.append(count)
            count = 0
            end_orf = 0
            star_orf = 1
        if aa == '*' and star_orf == 1 and end_orf == 0:
            if_orf.append(-1)
            num.append(count)
            count = 0
            end_orf = 1
            star_orf = 0
    num.append(count)

    score = cal_score(num, if_orf)

    return score


# 推荐算法
def cal_score(num, if_orf):
    return np.dot(np.array(num), np.array(if_orf))


# 推荐函数
def recommend(aas):
    r = {}
    for aa in aas:
        tempscore = num_tran_aa(aa)
        r.update({aa: tempscore})
    r_s = sorted(r.items(), key=lambda x: x[1], reverse=True)
    return r_s


def main():
    # Tkinter 文本框控件中第一个字符的位置是 1.0，可以用数字 1.0 或字符串"1.0"来表示。
    # "end"表示它将读取直到文本框的结尾的输入。我们也可以在这里使用 tk.END 代替字符串"end"
    dna = textExample.get('1.0', 'end')  # 获取文本框输入的内容
    dna = dna.upper().replace('U', 'T')
    dna = re.sub('[^ATGC]', '', dna)
    # dna = format_input()
    dna_t = dna[::-1].translate(str.maketrans(na))
    # print(dna)

    print('transl_table_1  Standard \n'
          'transl_table_2  Vertebrate mitochondrial \n'
          'transl_table_3  Yeast mitochondrial \n'
          'transl_table_4  Mold, protozoan and coelenterate mitochondrial, mycoplasma \n'
          'transl_table_5  lnvertebrate mitochondrial \n'
          'transl_table_6  Ciliate, dasycladacean and hexamita nuclear \n'
          'transl_table_7  Echinoderm and flatworm mitochondrial \n'
          'transl_table_8  Euplotid nuclear \n'
          'transl_table_9  Alternative yeast nuclear \n'
          'transl_table_10 Ascidian mitochondrial \n'
          'transl_table_11 Alternative flatworm mitochondrial \n'
          'transl_table_12 Chlorophycean mitochondrial \n'
          'transl_table_13 Trematode mitochondrial \n'
          'transl_table_14 Scenedesmus obliquus mitochondrial \n'
          'Choose transl_table number:')
    c = int(input())
    tran_aa_5to3 = tran(dna, transl_table[c - 1])
    tran_aa_3to5 = tran(dna_t, transl_table[c - 1])
    # for i in tran_aa_5to3 + tran_aa_3to5:
    #     print(i)

    result = recommend(tran_aa_5to3 + tran_aa_3to5)
    format_output(result)


if __name__ == '__main__':
    w = tk.Label(window, text='Please enter a protokaryotic DNA/RNA sequence'
                              ' - numbers and blanks are ignored:')
    w.pack()
    textExample = tk.Text(window, height=10)  # 文本输入框
    textExample.pack()  # 把Text放在window上面，显示Text这个控件

    # 按钮（#command绑定获取文本框内容的方法）
    btnRead = tk.Button(window, height=1, width=10, text='Translate', command=main)
    btnRead.pack()  # 显示按钮

    window.mainloop()  # 显示窗口
