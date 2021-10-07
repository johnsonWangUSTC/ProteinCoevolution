import Bio.PDB
import numpy as np


def calc_dist_matrix(model):

    residues = list()
    for r in model.get_residues():
        residues.append(r)
    L = len(residues)
    output = np.zeros([L, L], np.float)
    i = 0
    for r in residues:
        try:
            rb = r['CB']
        except KeyError:
            rb = r['CA']

        j = 0
        for s in residues:
            try:
                sb = s['CB']
            except KeyError:
                sb = s['CA']
            output[i, j] = sb - rb
            print(sb-rb, i, j)
            j += 1
        i += 1

    return output


def calc_contact_map(dist_matrix, gap=0):
    contact_map = np.zeros_like(dist_matrix)
    p = contact_map.shape[0]
    for i in range(p):
        for j in range(p):
            if 8 > dist_matrix[i, j] > 1e-8 and abs(i-j) > gap:
                contact_map[i, j] = 1 / (dist_matrix[i, j] ** 2)
            else:
                contact_map[i, j] = 0
    return contact_map


def merge_contact_map(a, b):

    ma = np.max(a)
    mb = np.max(b)
    return np.triu(a) + np.tril(b) * ma / mb


def to_contact(prec, gap=0):

    p = prec.shape[0]
    contact_map = np.zeros_like(prec)
    for i in range(p):
        for j in range(p):
            if abs(i-j) > gap and prec[i, j] != 0:
                contact_map[i, j] = np.abs(-prec[i, j] / prec[i, i] / prec[j, j])
            else:
                contact_map[i, j] = 0
    return contact_map


def top_k_match(contact_map, est_contact_map, k, gap=0):

    p = contact_map.shape[0]
    pos_1 = list()
    pos_2 = list()
    test_1 = np.triu(contact_map, gap+1)
    test_2 = np.triu(est_contact_map, gap+1)
    k = int(k)

    for i in range(k):
        p1 = np.argmax(test_1)
        r = p1 // p
        c = p1 - r * p
        test_1[r, c] = 0
        pos_1.append(p1)

        p2 = np.argmax(test_2)
        r = p2 // p
        c = p2 - r * p
        test_2[r, c] = 0
        pos_2.append(p2)

    return len(set(pos_1) & set(pos_2)) / k


def top_k_acc(contact_map, est_contact_map, k, gap=0):

    test = np.triu(est_contact_map, gap + 1)
    k = int(k)
    p = contact_map.shape[0]
    rec = 0

    for i in range(k):
        p1 = np.argmax(test)
        r = p1 // p
        c = p1 - r * p

        if test[r, c] == 0:
            return rec / i

        test[r, c] = 0
        if contact_map[r, c] != 0:
            rec += 1

    return rec / k


def apc(score):

    p = score.shape[0]
    row = np.mean(score,axis=1)
    m = np.mean(score)

    out = np.zeros_like(score)
    for i in range(p):
        for j in range(i+1):
            out[i, j] = out[j, i] = score[i, j] - row[i] * row[j] / m

    return out
