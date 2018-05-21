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
n = 5

H = ['h_'+str(i) for i in range(k)]
M = ['m_'+str(i) for i in range(m)]
N = ['n_'+str(i) for i in range(n)]

H_p = np.random.rand(k,2)
M_p = np.random.rand(m,2)
N_p = np.random.rand(n,2)

quantities = {f: 1 for i,f in enumerate(N)}
capacities = {h:(i%2 + 2)*2*2 for i,h in enumerate(H)}
type_dist = 'euclid'

vrp = VRPClass.VRP(H, N, H_p, N_p, quantities, capacities, type_dist, M = M, M_p = M_p)
#vrp.draw_problem()

h = 'h_0'
truck_capacity = capacities[h]
N
distance = vrp.distance_dict
values = lower_bound.possible_values(quantities,capacities[h])
values_pos = dict(zip(values,list(range(len(values)))))
quantities
direction = "left"

from cpp_lower_bounds import construct_q_paths
print(construct_q_paths(h,truck_capacity,N,distance,values,values_pos,quantities,direction))