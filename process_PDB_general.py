import Bio.PDB
import numpy as np
import pandas as pd
from sklearn.covariance import graphical_lasso
from plot_heatmap import *
from metrics import *

code = "1ctfA"
pdb_filename = "./data/pdb/"+code+".pdb"
gap = 23


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
    #return

    df = pd.read_csv("./data/corr/corr_"+code+".csv")
    corr = df.values[:, 1:]
    print(corr)
    heatmap(np.abs(corr), lim=[-1, 1])
    model = graphical_lasso(corr, alpha=.01, tol=1e-8, enet_tol=1e-8, max_iter=1000)
    prec = model[1]
    est_contact_map = to_contact(prec, gap=gap)
    #apc_score = np.abs(apc(est_contact_map))
    #heatmap(est_contact_map, lim='b', cmap='Blues')
    cpr = merge_contact_map(contact_map, est_contact_map)
    heatmap(cpr, lim='b', cmap='Blues')
    #print(cpr)
    par = [1, 2, 5, 10]
    for i in par:
        print("top-L/", i)
        print("new: ", top_k_acc(contact_map, est_contact_map, k=p/i, gap=gap))
        print("old: ", top_k_acc(contact_map, np.abs(corr), k=p/i, gap=gap))


if __name__ == '__main__':

    main()
