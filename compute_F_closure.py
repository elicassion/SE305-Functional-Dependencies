# Compute Functional Dependencies of F
# For each r in R find r+
# For each S in r+ output r->S

# Example:
# R = {A, B, C, D, E}
# F = {
#   A   ->  BC
#   CD  ->  E
#   B   ->  D
#   E   ->  A
# }

import os
import numpy as np
from itertools import product


def gen_subset(R):
    rsub = []
    n = len(R)
    if n == 0 or R == '\\emptyset':
        return ['']
    K = np.array(range(n))
    for i in range(2**n):
        e = list(bin(i))[2:]
        e = np.array(e) == '1'
        k = K[n-len(e):][e]
        g = [R[k_item] for k_item in k]
        rsub.append(g)
    return rsub

def attr_closure(r, F):
    res = [i for i in r]
    change = True
    while change:
        flag = False
        for k, v in F:
            if set(res) >= set(k):
                for i in v:
                    if i not in res:
                        res += i
                        flag = True
        change = flag
    return res

def compute_part_f_closure(r, F):
    rplus = attr_closure(r, F)
    s = gen_subset(rplus)
    return [(r, si) for si in s]

def get_superkeys(R, Fplus):
    superkeys = set([k for k in Fplus.keys() if R in Fplus[k]])
    return superkeys

def is_superkey(kb, R, Fplus):
    return R in Fplus[kb]

def is_candidate_key(k, R, Fplus):
    ksub = gen_subset(k)
    for kb in ksub:
        kb = norm_format(kb)
        if kb == k:
            continue
        if is_superkey(kb, R, Fplus):
            return False
    return True

def compute_candidate(R, Fplus):
    candidate_keys = []
    keys = get_superkeys(R, Fplus)
    # print (keys)
    for k in keys:
        if is_candidate_key(k, R, Fplus):
            candidate_keys.append(k)
    return candidate_keys


def norm_format(r):
    if r == [] or r == '':
        return '\\emptyset'
    else:
        q = list(r)
        q.sort()
        return ''.join(q)

def output(Fplus, candidate_keys, fc):
    fplusfile = open('Fclosure.txt', 'w')
    for k, v in Fplus:
        fplusfile.write("$%s \\rightarrow %s$\\\\\n" % (k, v))
    fplusfile.close()

    candidatefile = open('candidate_keys.txt', 'w')
    for k in candidate_keys:
        candidatefile.write("$%s$\\\\\n" % k)
    candidatefile.close()

    fcfile = open('canonical_cover.txt', 'w')
    for k,v in fc:
        fcfile.write("%s \\rightarrow %s, " % (k, v))
    fcfile.close()

def union_rule(Fc):
    fcdict = {}
    for k, v in Fc:
        if k not in fcdict.keys():
            fcdict[k] = []
        fcdict[k].append(v)
    nfcdict = {}
    for k, v in fcdict.items():
        nv = ''.join(set(''.join(v)))
        nfcdict[k] = nv
    return nfcdict.items()

def test_extraneous_attr_lhs(alpha, k, v, Fc):
    remain = ''.join(set(k) - set(alpha))
    # print ('remain:', remain)
    remain_closure = attr_closure(remain, Fc)
    if set(remain_closure) >= set(v):
        return True
    else:
        return False


def test_extraneous_attr_rhs(beta, k, v, Fc): 
    # print (Fc)
    # print (beta, k, v)
    Fdot = (set(Fc) - set([(k,v)])) | set([(k, ''.join(set(v)-set(beta)))])
    # print (Fdot)
    k_closure = attr_closure(k, Fdot)
    if beta in k_closure:
        return True
    else:
        return False


def delete_extraneous_attr(Fc):
    nfc = []
    change = False
    for k,v in Fc:
        nk = ''
        for alpha in k:
            if not test_extraneous_attr_lhs(alpha, k, v, Fc):
                nk += alpha
            else:
                change = True
        nv = ''
        for beta in v:
            if not test_extraneous_attr_rhs(beta, k, v, Fc):
                nv += beta
            else:
                change = True
        nfc.append((nk, nv))
    return nfc, change

def compute_canonical_cover(R, F):
    Fc = F
    change = True
    while change:
        Fc = union_rule(Fc)
        print (Fc)
        Fc, change = delete_extraneous_attr(Fc)
    return Fc


def main():
    R = 'ABCDE'
    F = [('A', 'BC'), 
         ('CD', 'E'),
         ('B',  'D'),
         ('E',  'A')]
    rsub = gen_subset(R)
    Fplus = set()
    for r in rsub:
        s = compute_part_f_closure(r, F)
        for k, v in s:
            k = norm_format(k)
            v = norm_format(v)
            Fplus.add((k, v))

    Fplus = list(Fplus)
    Fplus.sort(reverse = True)
    
    Fplusdict  = {}
    for k, v in Fplus:
        if k not in Fplusdict.keys():
            Fplusdict[k] = []
        Fplusdict[k].append(v)
    candidate_keys = compute_candidate(R, Fplusdict)

    Fc = compute_canonical_cover(R, F)
    output(Fplus, candidate_keys, Fc)

main()
