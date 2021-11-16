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
to_nums['EQUIV'] = to_nums['BLOSUM62']


def aa_to_num_McLachlan(aa):
    try:
        num = to_nums['McLachlan'][aa]
    except KeyError:
        num = -1
        #print("UNDEF\n")
    return num


def aa_to_num_BLOSUM62(aa):
    try:
        num = to_nums['BLOSUM62'][aa]
    except KeyError:
        num = -1
        print(aa, '\n')
    return num


def aa_to_num_EQUIV(aa):
    try:
        num = to_nums['EQUIV'][aa]
    except KeyError:
        num = -1
        print(aa, '\n')
    return num


def calc_corr(msa, type='McLachlan', gap_ratio=0.99):

    N = len(msa)
    L = len(list(msa[0]._seq))
    print(N, L)
    if type == 'McLachlan':
        aa_to_num = aa_to_num_McLachlan
        c_type = 0
    elif type == 'BLOSUM62':
        aa_to_num = aa_to_num_BLOSUM62
        c_type = 1
    elif type == 'EQUIV':
        aa_to_num = aa_to_num_EQUIV
        c_type = 2
    else:
        return
    count_gap = np.zeros(L)
    for i in range(N):
        seq = list(msa[i]._seq)
        for j in range(L):
            num = aa_to_num(seq[j])
            if num == -1:
                count_gap[j] += 1
    count_gap *= 1/N
    valid_pos = list()
    for j in range(L):
        if count_gap[j] <= gap_ratio:
            valid_pos.append(j)
    valid_L = len(valid_pos)
    print(count_gap, L, valid_L)

    T = N * valid_L
    MAT = c_int * T
    c_MSA = MAT()

    for i in range(N):
        seq = list(msa[i]._seq)
        for j in range(valid_L):
            num = aa_to_num(seq[j])
            c_MSA[i*valid_L+j] = num

    func = CDLL("calc_corr_C.dll")
    #func.calc_cov_C.argtypes = (MAT, c_int, c_int)
    func.calc_corr_C_v2.restype = (POINTER(POINTER(c_float)))
    c_corr = func.calc_corr_C_v2(c_MSA, c_int(N), c_int(valid_L), c_int(c_type))
    corr = np.zeros([L, L])
    for i in range(valid_L):
        for j in range(valid_L):
            corr[valid_pos[i], valid_pos[j]] = c_corr[i][j]
            print(c_corr[i][j])
    pop_idx = list()
    for i in range(valid_L):
        if corr[valid_pos[i], valid_pos[i]] == 0:
            pop_idx.append(i)
    pop_idx = sorted(pop_idx, reverse=True)
    print(pop_idx)
    for i in pop_idx:
        valid_pos.pop(i)
    corr = corr[valid_pos, :]
    corr = corr[:, valid_pos]
    return corr, valid_pos


def main():

    code = "1ag6A"
    typ = "EQUIV"
    msa = Bio.AlignIO.read("data/msa/"+code+".fasta", "fasta")

    corr, valid_pos = calc_corr(msa, type=typ, gap_ratio=0.1)
    print("Calculation Completed\n")
    np.savetxt("data/corr/"+typ+"_v2/corr_"+typ+"_"+code+".txt", corr, delimiter=" ")
    np.savetxt("data/corr/" + typ + "_v2/vpos_" + typ + "_" + code + ".txt", valid_pos, delimiter=" ")
    heatmap(corr, lim='a')


if __name__ == '__main__':
    main()