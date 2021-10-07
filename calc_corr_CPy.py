import Bio.AlignIO
import numpy as np
import pandas as pd
from ctypes import *
from plot_heatmap import *

to_names = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",  "-"]
to_nums = dict()
to_nums['McLachlan'] = {'A':0, 'C':1, 'D':2, 'E':3, 'F':4, 'G':5, 'H':6, 'I':7, 'K':8, 'L':9, 'M':10, 'N':11, 'P':12,
           'Q':13, 'R':14, 'S':15, 'T':16, 'V':17, 'W':18, 'Y':19}
to_nums['BLOSUM62'] = {'A':0, 'R':1, 'N':2, 'D':3,  'C':4,  'Q':5,  'E':6,  'G':7,  'H':8,  'I':9,  'L':10,  'K':11,  'M':12,
                       'F':13, 'P':14,  'S':15,  'T':16,  'W':17,  'Y':18,  'V':19,  'B':20,  'Z':21,  'X':22, '-':23}


def aa_to_num_McLachlan(aa):
    try:
        num = to_nums['McLachlan'][aa]
    except KeyError:
        num = -1
    return num


def aa_to_num_BLOSUM62(aa):
    try:
        num = to_nums['BLOSUM62'][aa]
    except KeyError:
        num = -1
    return num




mat_mclachlan = np.array([[8, 1, 3, 4, 1, 3, 3, 2, 3, 2, 3, 3, 4, 3, 2, 4, 3, 3, 1, 1],
                          [1, 9, 1, 0, 0, 1, 3, 1, 0, 0, 3, 1, 0, 0, 1, 2, 2, 1, 2, 1],
                          [3, 1, 8, 5, 1, 3, 4, 0, 3, 1, 2, 5, 3, 4, 1, 3, 3, 1, 0, 1],
                          [4, 0, 5, 8, 0, 3, 2, 1, 4, 1, 1, 4, 4, 5, 3, 4, 4, 2, 1, 2],
                          [1, 0, 1, 0, 9, 0, 4, 3, 0, 5, 5, 0, 1, 0, 1, 2, 1, 3, 6, 6],
                          [3, 1, 3, 3, 0, 8, 2, 1, 3, 1, 1, 3, 3, 2, 3, 3, 2, 2, 1, 0],
                          [3, 3, 4, 2, 4, 2, 8, 2, 4, 2, 3, 4, 3, 4, 5, 3, 4, 2, 3, 4],
                          [2, 1, 0, 1, 3, 1, 2, 8, 1, 5, 5, 1, 1, 0, 1, 2, 3, 5, 3, 3],
                          [3, 0, 3, 4, 0, 3, 4, 1, 8, 2, 1, 4, 3, 4, 5, 3, 3, 2, 1, 1],
                          [2, 0, 1, 1, 5, 1, 2, 5, 2, 8, 6, 1, 1, 3, 2, 2, 3, 5, 3, 3],
                          [3, 3, 2, 1, 5, 1, 3, 5, 1, 6, 8, 2, 1, 3, 1, 2, 3, 4, 1, 2],
                          [3, 1, 5, 4, 0, 3, 4, 1, 4, 1, 2, 8, 1, 4, 3, 5, 3, 1, 0, 2],
                          [4, 0, 3, 4, 1, 3, 3, 1, 3, 1, 1, 1, 8, 3, 3, 3, 3, 2, 0, 0],
                          [3, 0, 4, 5, 0, 2, 4, 0, 4, 3, 3, 4, 3, 8, 5, 4, 3, 2, 2, 1],
                          [2, 1, 1, 3, 1, 3, 5, 1, 5, 2, 1, 3, 3, 5, 8, 4, 3, 2, 3, 2],
                          [4, 2, 3, 4, 2, 3, 3, 2, 3, 2, 2, 5, 3, 4, 4, 8, 5, 2, 3, 3],
                          [3, 2, 3, 4, 1, 2, 4, 3, 3, 3, 3, 3, 3, 3, 3, 5, 8, 3, 2, 1],
                          [3, 1, 1, 2, 3, 2, 2, 5, 2, 5, 4, 1, 2, 2, 2, 2, 3, 8, 2, 3],
                          [1, 2, 0, 1, 6, 1, 3, 3, 1, 3, 1, 0, 0, 2, 3, 3, 2, 2, 9, 6],
                          [1, 1, 1, 2, 6, 0, 4, 3, 1, 3, 2, 2, 0, 1, 2, 3, 1, 3, 6, 9]])


def calc_corr(msa, type='McLachlan'):

    N = len(msa)
    L = len(list(msa[0]._seq))
    print(N, L)
    T = N * L
    MAT = c_int * T
    c_MSA = MAT()
    corr = np.zeros([L, L])
    if type == 'McLachlan':
        aa_to_num = aa_to_num_McLachlan
        c_type = 0
    elif type == 'BLOSUM62':
        aa_to_num = aa_to_num_BLOSUM62
        c_type = 1
    else:
        return
    for i in range(N):
        seq = list(msa[i]._seq)
        for j in range(L):
            c_MSA[i*L+j] = aa_to_num(seq[j])
    func = CDLL("calc_corr_C.dll")
    #func.calc_corr_C.argtypes = (MAT, c_int, c_int)
    func.calc_corr_C.restype = (POINTER(POINTER(c_float)))
    c_corr = func.calc_corr_C(c_MSA, c_int(N), c_int(L), c_int(c_type))
    corr = np.zeros([L, L])
    for i in range(L):
        for j in range(L):
            corr[i, j] = c_corr[i][j]
            print(c_corr[i][j])
    return corr



def main():

    code = "1g2rA"
    typ = "BLOSUM62"
    msa = Bio.AlignIO.read("data/msa/"+code+".fasta", "fasta")
    '''
    COL = c_double * 3
    a = COL
    b = COL
    a[0] = 1
    '''
    corr = calc_corr(msa, type="BLOSUM62")
    np.savetxt("data/corr/"+typ+"/corr_"+typ+"_"+code+".txt", corr, delimiter=" ")
    heatmap(corr, lim='a')


if __name__ == '__main__':
    main()
