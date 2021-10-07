import Bio.PDB
import numpy as np
from plot_heatmap import *

pdb_code = "1guu"
pdb_filename = "1guuA.pdb"


def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector = residue_one["CA"].coord - residue_two["CA"].coord
    return np.sqrt(np.sum(diff_vector * diff_vector))


def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = np.ones((len(chain_one), len(chain_two)), np.float) * -1
    for row, residue_one in enumerate(chain_one):
        for col, residue_two in enumerate(chain_two):

            if residue_one != residue_two:
                # compute distance between CA atoms
                try:
                    distance = residue_one['CB'] - residue_two['CB']
                except KeyError:
                    ## no CA atom, e.g. for H_NAG
                    continue
                if distance < 8:
                    print(residue_one, residue_two, distance)
                answer[row, col] = distance
    return answer


def main():

    structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
    print(structure)
    model = structure[0]
    chain = model['A']
    print(chain)

    dist_matrix = calc_dist_matrix(model['A'], model['A'])
    dist_matrix = dist_matrix[0:50, 0:50]
    p = len(dist_matrix[0, ])
    print(p)
    contact_map = np.zeros_like(dist_matrix)
    for i in range(p):
        for j in range(p):
            if 8 > dist_matrix[i, j] > 0 and abs(i-j) > 8:
                contact_map[i, j] = 1
    heatmap(contact_map)
    #heatmap(dist_matrix)


if __name__ == '__main__':

    main()
