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


from cpp_lower_bounds import lower_bound as lb
#lower_bound.lower_bound()

lamb = {}
mu = {}
for j in N:
    lamb[j] = 0
for h in H:
    mu[h] = 0
#print(lower_bound.lower_bound(H,capacities,N,quantities,distance,mu,lamb)[0])

H_ = (np.array(range(len(H)))+len(N)).astype(int)
N_ = np.array(range(len(N))).astype(int)
capacities_ = np.array([capacities[h] for h in H]).astype(int)
quantities_ = np.array([quantities[n] for n in N]).astype(int)
#print(vrp.distance_dict)
distance_ = vrp.distance.astype('float64')
#print(distance_)

mu_ = np.array([mu[h] for h in H]).astype('float64')
lamb_ = np.array([lamb[n] for n in N]).astype('float64')
#print(len(H_))
#print(len(N_))

result = lb(H_,capacities_,N_,quantities_,distance_,mu_,lamb_)
print (result['z_lb'])
z_lb, theta, rho, u = lower_bound.lower_bound(H,capacities,N,quantities,vrp.distance_dict,mu,lamb)
print(z_lb)
