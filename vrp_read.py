import sys
import VRPClass_light as VRPClass

import cpp_lower_bounds
import time
import os


import numpy as np
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist
import copy
from copy import deepcopy
import lower_bound
lower_bound = reload(lower_bound)

import pickle

time_limit = 120

start = time.time()
stats = pickle.load(open( "stats.p",'rb' ) )
end = time.time()

filename = 'results.txt'

k = int(stats["size"]["k"])
n = int(stats["size"]["n"])
m = int(stats["size"]["m"])
capacity = stats["problem"]["capacity"]

H = ['h_'+str(i) for i in range(k)]
N = ['n_'+str(i) for i in range(n)]
M = ['m_'+str(i) for i in range(m)]

quantities = {f: 1 for i,f in enumerate(N)}
capacities = {h:capacity for i,h in enumerate(H)}

type_dist = 'euclid'
vrp = VRPClass.VRP(H, N, stats["problem"]["H_p"], stats["problem"]["N_p"], quantities, capacities, type_dist, M = M, M_p = stats["problem"]["M_p"])


for i in range(len(stats["result"])):
    routes = lower_bound.primal_solver(stats["result"][i],len(N),H, quantities, capacities, time_limit)
    ub = vrp.total_distance(routes)
    lb = stats["result"][i]["z_lb"]
    gap = (ub-lb)/ub
    gamma_guarantee = stats["result"][i]["gamma_guarantee"]
    if os.path.exists(filename):
        append_write = 'a' # append if already exists
    else:
        append_write = 'w' # make a new file if not
    with open(filename, append_write) as file:
        file.write(
            str(stats["size"]["n"]) +
            ";" + str(stats["size"]["k"]) +
            ";" + str(stats["problem"]["capacity"]) +
            ";" + str(stats["time"]/(len(stats["result"])-1)) +
            ";" + str(i) +
            ";" + str(ub) +
            ";" + str(lb) +
            ";" + str(gap) +
            ";" + str(gamma_guarantee) + '\n')
