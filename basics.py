import os,sys
# from math import *
import math

# I interpret bond lengths below bond_dist2 as "bonded"
bond_dist1 = 1.5
bond_dist2 = 1.76
### subject to refinement!
too_close_dist = 0.99
# too_close_dist = 0.74 (number printed in R160 runs... Gaussian specific)

# How far out from structure to put the oxygen.
ox_height = 1.25;  # Oxygen places at midpoint of 5-6 bond, 
# expanded outward from center of structure by this factor
### subject to  ref...!

indene_height = 1.54
indene_carbon_height = 1.1

# How much outward carbons get pushed when oxygen arrives
carbon_expansion_factor = 1.4


def vec_add(x,y):
    return [x[i] + y[i] for i in range(0,len(x))]

def vec_sub(x,y):
    return [x[i] - y[i] for i in range(0,len(x))]

def vec_ave(x,y):
    return [0.5 * (x[i] + y[i]) for i in range(0,len(x))]

def vec_axpy(a,x,y):
    return [a * x[i] + y[i] for i in range(0,len(x))]

def vec_norm(x):
    return math.sqrt(sum([x[i]*x[i] for i in range(0,len(x))]))

def vec_wtnorm(w,x):
    return math.sqrt(sum([w[i]*w[i]*x[i]*x[i] for i in range(0,len(x))]))

def vec_cross(x,y):
    return [x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0]]

def vec_dot(x,y):
    return sum([x[i]*y[i] for i in range(0,len(x))]) 

def vec_cos(x,y):
    return (vec_dot(x,y)/(vec_norm(x) * vec_norm(y)))

def vec_normalize(x):
    d = vec_norm(x)
    return ([x[i]/d for i in range(0,len(x))])

def vec_normal_to(p):
    v1 = vec_sub(p[1], p[0])
    v2 = vec_sub(p[2], p[0])
    n1 = vec_cross(v1,v2)
    return n1

def coplanar(xyz,w):
    """
    Given A,B,C,D, check that (B-A) X (C-A) * (D-A) ~ 0.
    lack of _exact_ planarity in planes of interest requires a "slush" parameter
    """
    v1 = vec_sub(xyz[1], xyz[0])
    v2 = vec_sub(xyz[2], xyz[0])
    v3 = vec_sub(w, xyz[0])
    n1 = vec_cross(v1,v2)
    d1 = vec_dot(v3, n1)
    # this parameter has to be tuned until you get the right answer:!
    eps=0.6
    if (abs(d1) < eps):
        return True
    else:
        return False

def is_subset(test_idx, idx):
    subset = True
    for atm_idx in test_idx:
        if (not atm_idx in idx):
            subset = False
    return (subset)


def get_subset(test_idx, idx):
    """
    already know test_idx is subset of idx, not need indices in idx that match
    sort of a hack right now to store this info.
    """
    idx_inter = []
    for i in range(0,len(idx)):
        if idx[i] in test_idx:
            idx_inter.append(i)
    return (idx_inter)


def mat_vec_mult(M,v):
    z = [0 for i in range(0,len(M))]
    for i in range(0,len(z)):
        for j in range(0,len(v)):
            z[i] += M[i][j] * v[j]
    return z

def mat_zero(n,m):
    A = [[0 for i in range(m)] for j in range(n)]
    return A

def matT_mat_mult(A,B):
    # compute A^T B
    n = len(A)
    m = len(A[0])
    p = len(B[0])
    assert(len(B) == n)
    C = mat_zero(m,p)
    for i in range(m):
        for j in range(p):
            for k in range(n):
                C[i][j] += A[k][i] * B[k][j]
    return C

def mat_mat_mult(A,B):
    # compute A B
    n = len(A)
    m = len(A[0])
    p = len(B[0])
    assert(len(B) == m)
    C = mat_zero(n,p)
    for i in range(n):
        for j in range(p):
            for k in range(m):
                C[i][j] += A[i][k] * B[k][j]
    return C


def mat3_det(A):
    a = A[0][0] * (A[1][1]*A[2][2]-A[1][2]*A[2][1])
    b = A[0][1] * (A[1][0]*A[2][2]-A[1][2]*A[2][0])
    c = A[0][2] * (A[1][0]*A[2][1]-A[2][0]*A[1][1])
    return (a-b+c)

def mat3_trace(A):
    return (A[0][0] + A[1][1] + A[2][2])

