import Bio.PDB
import numpy as np
import pandas as pd
from sklearn.covariance import graphical_lasso
from plot_heatmap import *
from metrics import *

pdb_code = "1brfA"
pdb_filename = "./data/pdb/1ag6A.pdb"
gap = 8


def main():

    structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
    print(structure)
    model = structure[0]
    print(type(model))
    print(dir(model))
    print(vars(model))

    dist_matrix = calc_dist_matrix(model)
    #print(dist_matrix)
    p = len(dist_matrix[0, ])
    print(p)
    contact_map = calc_contact_map(dist_matrix, gap=gap)
    #heatmap(contact_map, lim='b', cmap='Blues')
    #heatmap(dist_matrix)


    df = pd.read_csv("./data/corr/corr_1ag6A.csv")
    corr = df.values[:, 1:]
    #heatmap(np.abs(corr))
    model = graphical_lasso(corr, alpha=.1, tol=1e-8, enet_tol=1e-5, max_iter=500)
    prec = model[1]
    est_contact_map = to_contact(prec, gap=gap)
    #heatmap(est_contact_map, lim='b', cmap='Blues')
    cpr = merge_contact_map(contact_map, est_contact_map)
    #heatmap(cpr, lim='b', cmap='Blues')
    #print(cpr)
    print(top_k_acc(contact_map, est_contact_map, k=p, gap=gap))
    print(top_k_acc(contact_map, np.abs(corr), k=p, gap=gap))


if __name__ == '__main__':

    main()
