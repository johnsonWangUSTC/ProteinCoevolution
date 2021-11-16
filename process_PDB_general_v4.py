import Bio.PDB
import numpy as np
import pandas as pd
from sklearn.covariance import graphical_lasso
from plot_heatmap import *
from metrics import *


code = "1dlwA"
typ = "EQUIV"
pdb_filename = "./data/pdb/"+code+".pdb"
gap = 4


def main():

    structure = Bio.PDB.PDBParser().get_structure(code, pdb_filename)
    print(structure)
    model = structure[0]

    dist_matrix = calc_dist_matrix(model)
    #print(dist_matrix)
    p = len(dist_matrix[0, ])
    print(p)
    contact_map = calc_contact_map(dist_matrix, gap=gap)

    #heatmap(contact_map, lim='b', cmap='Blues')
    #heatmap(dist_matrix)

    cov = np.loadtxt("data/cov/"+typ+"/cov_"+typ+"_"+code+".txt")
    valid_pos = np.loadtxt("data/cov/" + typ + "/vpos_" + typ + "_" + code + ".txt")
    valid_pos = np.array(valid_pos, dtype=np.int)
    print(valid_pos)

    corr = np.loadtxt("data/corr/" + typ + "/corr_" + typ + "_" + code + ".txt")
    #heatmap(corr)
    model = graphical_lasso(cov, alpha=1e-4, tol=1e-8, enet_tol=1e-8, max_iter=500)
    prec = model[1]

    ext_prec = np.eye(p)
    valid_L = len(valid_pos)
    for i in range(valid_L):
        for j in range(valid_L):
            ext_prec[valid_pos[i], valid_pos[j]] = prec[i, j]
    print(np.max(ext_prec))
    print(np.linalg.eigvals(prec))
    est_contact_map = to_contact(ext_prec, gap=gap)

    cpr = merge_contact_map(est_contact_map, contact_map)
    heatmap(cpr, lim='b', cmap='Blues')

    par = [1, 2, 5, 10]
    for i in par:
        print("top-L/", i)
        print("new: ", top_k_acc(contact_map, est_contact_map, k=p / i, gap=gap))
        print("old: ", top_k_acc(contact_map, np.abs(corr), k=p / i, gap=gap))
        print("more: ", top_k_acc(contact_map, to_contact(np.linalg.inv(corr), gap=gap), k=p / i, gap=gap))


if __name__ == '__main__':

    main()
