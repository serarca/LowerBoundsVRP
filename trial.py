# In this code we implement the lower bounds of Baldacci's implementation
# Import VRP package
import sys
sys.path.insert(0, '/Users/sergiocamelo/Dropbox/Sergio-Joann/Code')
import VRPClass

import numpy as np
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist
import copy
from copy import deepcopy
import lower_bound
lower_bound = reload(lower_bound)

# We want to solve the Vehicle Routing Problem
# We first generate a random problem
k = 2
m = 1
n = 10

H = ['h_'+str(i) for i in range(k)]
M = ['m_'+str(i) for i in range(m)]
N = ['n_'+str(i) for i in range(n)]

H_p = np.random.rand(k,2)
M_p = np.random.rand(m,2)
N_p = np.random.rand(n,2)

quantities = {f: 1 for i,f in enumerate(N)}
capacities = {h:5 for i,h in enumerate(H)}
type_dist = 'euclid'

vrp = VRPClass.VRP(H, N, H_p, N_p, quantities, capacities, type_dist, M = M, M_p = M_p)
#vrp.draw_problem()


import cpp_lower_bounds
#lower_bound.lower_bound()


#lower_bound.optimize_lower_bound_c(100, 10, 0.1, H,capacities,N,quantities,vrp.distance)

#lower_bound.optimize_lower_bound(100, 10, 0.1, H,capacities,N,quantities,vrp.distance_dict)

Delta = 10;
gamma = 10;
h = 'h_0'
h_ = len(N);
capacity = capacities[h]
distance_mat = vrp.distance;
direction = "left"

path_c = lower_bound.GENPATH_c(Delta, gamma, h_, capacity, N, quantities, distance_mat, direction)
path = lower_bound.GENPATH(Delta, gamma, h, capacity, N, quantities, vrp.distance_dict, direction)

print(path_c)
print(" ")
#print(path)
#print(path == path_c)

route_c = lower_bound.GENROUTE_c(Delta, gamma, h_, capacity, N, quantities, distance_mat, distance_mat)
route = lower_bound.GENROUTE(Delta, gamma, h, capacity, N, quantities, vrp.distance_dict)

print("route_c")
print(route_c)
print(" ")
print("route")
print(route)
print(route == route_c)
