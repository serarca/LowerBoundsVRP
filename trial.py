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
capacities = {h:5 for i,h in enumerate(H)}
type_dist = 'euclid'

vrp = VRPClass.VRP(H, N, H_p, N_p, quantities, capacities, type_dist, M = M, M_p = M_p)
#vrp.draw_problem()


from cpp_lower_bounds import lower_bound as lb
#lower_bound.lower_bound()

lamb = {}
mu = {}
for j in N:
    lamb[j] = 0
for h in H:
    mu[h] = 0
#print(lower_bound.lower_bound(H,capacities,N,quantities,distance,mu,lamb)[0])

H_ = np.array(range(len(H)))
N_ = np.array(range(len(N)))
capacities_ = [capacities[h] for h in H]
quantities_ = [quantities[n] for n in N]
print(vrp.distance_dict)
distance_ = vrp.distance
print(distance_)

mu_ = [mu[h] for h in H]
lamb_ = [lamb[n] for n in N]
print(len(H_))
print(len(N_))

result = lb(H_,capacities_,N_,quantities_,distance_,mu_,lamb_)
print (result.z_lb)
