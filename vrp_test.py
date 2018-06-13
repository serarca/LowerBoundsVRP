import sys
import VRPClass_light as VRPClass

import cpp_lower_bounds
import time


import numpy as np
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist
import copy
from copy import deepcopy
import lower_bound
lower_bound = reload(lower_bound)

# We want to solve the Vehicle Routing Problem
# We first generate a random problem
k = int(sys.argv[1])
n = int(sys.argv[2])
capacity = int(sys.argv[3])
z_ub = float(sys.argv[4])


Delta = 2000
gamma = 1
iterations_m2 = 5




# Parameters for the search algorithm that will not change
m = 1
iterations_grad_m1 = 200
iterations_grad_m2 = 100
Delta_zero = 1000
Delta_final = Delta
gamma_zero = -10**(-14)
gamma_final = gamma
epsilon = 0.1



## Generate random coord
H = ['h_'+str(i) for i in range(k)]
M = ['m_'+str(i) for i in range(m)]
N = ['n_'+str(i) for i in range(n)]

H_p = np.random.rand(k,2)
M_p = np.random.rand(m,2)
N_p = np.random.rand(n,2)

quantities = {f: 1 for i,f in enumerate(N)}
capacities = {h:capacity for i,h in enumerate(H)}
type_dist = 'euclid'
vrp = VRPClass.VRP(H, N, H_p, N_p, quantities, capacities, type_dist, M = M, M_p = M_p)



geo_distance = vrp.distance


H_ = (np.array(range(len(H)))+len(N)).astype(int)
N_ = np.array(range(len(N))).astype(int)
capacities_ = np.array([capacities[h] for h in H]).astype(int)
quantities_ = np.array([quantities[n] for n in N]).astype(int)
geo_distance_ = geo_distance.astype("float64")


start = time.time()
# We time the problem
result = lower_bound.construct_lower_bound_c(iterations_grad_m1,iterations_grad_m2,iterations_m2,z_ub,Delta,Delta_zero,Delta_final,gamma,gamma_zero,gamma_final,epsilon,H,capacities,N,quantities,geo_distance)
end = time.time()
total_time = end - start


stats = {
    "time":total_time,
    "size":{"k" : len(H), "m" : len(M), "n" : len(N)},
    "solver":{
        "iterations_grad_m1" : iterations_grad_m1,
        "iterations_grad_m2" : iterations_grad_m2,
        "iterations_m2" : iterations_m2,
        "z_ub" : z_ub,
        "Delta" : Delta,
        "Delta_zero" : Delta_zero,
        "Delta_final" : Delta_final,
        "gamma" : gamma,
        "gamma_zero" : gamma_zero,
        "gamma_final" : gamma_final,
        "epsilon" : epsilon},
    "result":result,
    "problem":{"H_p":H_p,
        "M_p":M_p,
        "N_p":N_p,
        "capacity":capacity}
    }
import pickle
pickle.dump(stats, open( "stats.p", "wb" ) )
