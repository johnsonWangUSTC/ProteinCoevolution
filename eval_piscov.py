import Bio.PDB
import numpy as np
import pandas as pd
from sklearn.covariance import graphical_lasso
from plot_heatmap import *
from metrics import *
from ctypes import *


code = "1ctfA"
pdb_filename = "./data/pdb/"+code+".pdb"
gap = 23


def main():

    with open("data/psicov/" + code + ".txt", "r") as f:
        txt = f.readlines()[0]

    data = txt.split(sep='\t')
    print(data)
    M = int(((len(data))-1) / 2)
    print(M)
    idx = np.zeros([2, M], dtype=np.int)
    for i in range(M):
        idx[0, i] = int(data[2 * i])
        idx[1, i] = int(data[2 * i + 1])
    print(idx)

    structure = Bio.PDB.PDBParser().get_structure(code, pdb_filename)
    model = structure[0]

    dist_matrix = calc_dist_matrix(model)
    #print(dist_matrix)
    p = len(dist_matrix[0, ])
    print(p)
    contact_map = calc_contact_map(dist_matrix, gap)


    par = [1, 2, 5, 10]
    for i in par:
        print("top-L/", i)
        print("PSICOV: ", top_k_acc_psicov(contact_map, idx, k=p / i, gap=gap))


if __name__ == '__main__':

    main()
