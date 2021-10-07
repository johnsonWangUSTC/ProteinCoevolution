import Bio.AlignIO
import numpy as np
import pandas as pd

to_names = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",  "-"]
to_nums = {'A':0, 'C':1, 'D':2, 'E':3, 'F':4, 'G':5, 'H':6, 'I':7, 'K':8, 'L':9, 'M':10, 'N':11, 'P':12,
           'Q':13, 'R':14, 'S':15, 'T':16, 'V':17, 'W':18, 'Y':19,  '-':20}


def seq_to_nums(seq):

    L = len(seq)
    nums = np.zeros(L, dtype=np.int)
    for i in range(L):
        try:
            nums[i] = to_nums[seq[i]]
            #print(nums[i])
        except KeyError:
            nums[i] = 20
            #print(20)
    return nums


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


def calc_corr(msa, gap_ratio=0.99):

    if (gap_ratio < 0) or (gap_ratio > 1):

        print("gap_ratio must be in the [0,1] range.")
        return
    N = len(msa)
    L = len(list(msa[0]._seq))
    print(N, L)
    MSA = np.zeros([N, L])
    corr = np.zeros([L, L])
    for i in range(N):
        seq = list(msa[i]._seq)
        MSA[i, ] = np.array(seq_to_nums(seq),dtype=np.int)
    print(MSA)

    dat = np.zeros([N*N, L])
    for k in range(L):
        print(k)
        for i in range(N):
            for j in range(i+1):
                ai = int(MSA[i, k])
                aj = int(MSA[j, k])
                print()
                if ai != 20 and aj != 20:
                    dat[(i-1)*N + j, k] = dat[(j-1)*N + i, k] = mat_mclachlan[ai, aj]
                else:
                    dat[(i-1)*N + j, k] = dat[(j-1)*N + i, k] = 255

    for i in range(L):
        for j in range(i+1):
            print(i, j)
            bi = dat[:, i]
            pi = np.where(bi != 255)
            bj = dat[:, j]
            pj = np.where(bj != 255)
            pos = np.intersect1d(pi, pj)
            bi = bi[pos]
            bj = bj[pos]
            corr[i, j] = corr[j, i] = np.corrcoef(bi, bj)[0, 1]

    return corr


def main():

    msa = Bio.AlignIO.read("1ag6A.fasta", "fasta")
    msa = msa[0:100]
    corr = calc_corr(msa)
    print(corr)
    np.savetxt("corr.txt", corr, delimiter=" ")
    corr = np.loadtxt("corr.txt")

    print(corr)




if __name__ == '__main__':
    main()
