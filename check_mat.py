import numpy as np
from plot_heatmap import *
from metrics import *
import Bio.AlignIO

typ = "EQUIV"
code = "1aapA"
gap = 4
pdb_filename = "./data/pdb/"+code+".pdb"


def main():
    cov = np.loadtxt("data/cov/" + typ + "/cov_" + typ + "_" + code + ".txt")
    valid_pos = np.loadtxt("data/cov/" + typ + "/vpos_" + typ + "_" + code + ".txt")
    valid_pos = np.array(valid_pos, dtype=np.int)

    structure = Bio.PDB.PDBParser().get_structure(code, pdb_filename)
    model = structure[0]
    msa = Bio.AlignIO.read("data/msa/" + code + ".fasta", "fasta")

    N = len(msa)
    L = len(list(msa[0]._seq))

    ext_cov = np.zeros([L, L])
    valid_L = len(valid_pos)
    for i in range(valid_L):
        for j in range(valid_L):
            ext_cov[valid_pos[i], valid_pos[j]] = cov[i, j]

    dist_matrix = calc_dist_matrix(model)
    contact_map = calc_contact_map(dist_matrix, gap=gap)
    cpr = merge_contact_map(np.abs(ext_cov), contact_map)
    heatmap(cpr, lim='b', cmap='Blues')

    # heatmap(contact_map, lim='b', cmap='Blues')
    # heatmap(dist_matrix)

if __name__ == '__main__':
    main()