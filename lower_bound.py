# In this code we implement the lower bounds of Baldacci's implementation
# Import VRP package
import sys
sys.path.insert(0, '/Users/sergiocamelo/Dropbox/Sergio-Joann/Code')
import VRPClass
import ipdb

import numpy as np
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist
import copy
from copy import deepcopy
import cpp_lower_bounds



def possible_values(quantities, maximum):
    possible_values = []
    q_int = [int(quantities[n]) for n in quantities.keys()]
    def GCD(a, b):
        if b == 0:
            return a
        else:
            return GCD(b, a % b)
    reducer = int(reduce(GCD, (q_int)))
    s = reducer
    while s<=maximum:
        possible_values.append(s)
        s+=reducer
    return possible_values

# Generates a dictionary of reduced costs given a distance and a lambda vector
def reduced_cost_dict(lamb, distance_dictionary, N):
    distance = copy.deepcopy(distance_dictionary)
    for k1 in distance.keys():
        for k2 in distance[k1].keys():
            if k1 in N:
                distance[k1][k2] = distance[k1][k2] - 1.0/2*lamb[k1]
            if k2 in N:
                distance[k1][k2] = distance[k1][k2] - 1.0/2*lamb[k2]
    return distance

# Generates a matrix of reduced costs given a distance and a lambda vector
def reduced_cost_mat(lamb_, distance_mat, N_):
    distance = copy.deepcopy(distance_mat)
    for i in range(len(N_)):
        for j in range(distance_mat.shape[0]):
            distance[i][j] -= 1.0/2*lamb_[i]
    for i in range(len(N_)):
        for j in range(distance_mat.shape[0]):
            distance[j][i] -= 1.0/2*lamb_[i]
    return distance.astype('float64')

def construct_q_paths(h,truck_capacity,N,distance,values,values_pos,quantities,direction):

    # Initialize the routes
    f = {}; phi = {}; p = {}; q_route = {}; q_route_2 = {};
    for l in range(len(values)):
        f[l] = {}; phi[l] = {}; p[l] = {}; q_route[l] = {}; q_route_2[l] = {};
        for n in N:
            f[l][n] = float('inf'); phi[l][n] = float('inf'); p[l][n] = float('inf');
            q_route[l][n] = []; q_route_2[l][n] = [];

    # Initialize the routes
    for n in N:
        q = quantities[n]
        if q<= truck_capacity:
            l = values_pos[q]
            if direction == 'left':
                f[l][n] = distance[h][n]
            else:
                f[l][n] = distance[n][h]
            p[l][n] = h
            q_route[l][n] = [h, n]

    # Calculate the recursion
    for l in range(len(values)):
        Q = values[l]
        g = {}
        g_type = {}
        for x_i in N:
            q_p = Q - quantities[x_i]
            g[x_i] = {}
            g_type[x_i] = {}
            if q_p>0:
                l_p = values_pos[q_p]
                for x_j in N:
                    if x_i!=x_j:
                        if p[l_p][x_j]!=x_i:
                            g[x_i][x_j] = f[l_p][x_j] + distance[x_j][x_i]
                            g_type[x_i][x_j] = q_route[l_p][x_j] + [x_i]
                        else:
                            g[x_i][x_j] = phi[l_p][x_j] + distance[x_j][x_i]
                            g_type[x_i][x_j] = q_route_2[l_p][x_j] + [x_i]
        for x_i in N:
            q_p = Q - quantities[x_i]
            if q_p > 0:
                arg_min_1 = min(g[x_i], key=g[x_i].get)
                p[l][x_i] = arg_min_1
                f[l][x_i] = g[x_i][arg_min_1]
                del g[x_i][arg_min_1]
                arg_min_2 = min(g[x_i], key=g[x_i].get)
                phi[l][x_i] = g[x_i][arg_min_2]
                q_route[l][x_i] = g_type[x_i][arg_min_1]
                q_route_2[l][x_i] = g_type[x_i][arg_min_2]
    return (f, phi, p, q_route, q_route_2)

def construct_q_paths_debug(h,truck_capacity,N,distance,values,values_pos,quantities,direction):

    # Initialize the routes
    f = {}; phi = {}; p = {}; q_route = {}; q_route_2 = {};
    for l in range(len(values)):
        f[l] = {}; phi[l] = {}; p[l] = {}; q_route[l] = {}; q_route_2[l] = {};
        for n in N:
            f[l][n] = float('inf'); phi[l][n] = float('inf'); p[l][n] = float('inf');
            q_route[l][n] = []; q_route_2[l][n] = [];

    # Initialize the routes
    for n in N:
        q = quantities[n]
        if q<= truck_capacity:
            l = values_pos[q]
            if direction == 'left':
                f[l][n] = distance[h][n]
            else:
                f[l][n] = distance[n][h]
            p[l][n] = h
            q_route[l][n] = [h, n]

    # Calculate the recursion
    for l in range(len(values)):
        Q = values[l]
        g = {}
        g_type = {}
        for x_i in N:
            q_p = Q - quantities[x_i]
            g[x_i] = {}
            g_type[x_i] = {}
            if q_p>0:
                l_p = values_pos[q_p]
                for x_j in N:
                    if x_i!=x_j:
                        if p[l_p][x_j]!=x_i:
                            g[x_i][x_j] = f[l_p][x_j] + distance[x_j][x_i]
                            g_type[x_i][x_j] = (0,l_p)
                        else:
                            g[x_i][x_j] = phi[l_p][x_j] + distance[x_j][x_i]
                            g_type[x_i][x_j] = (1,l_p)
        for x_i in N:
            q_p = Q - quantities[x_i]
            if q_p > 0:
                min_1 = float('inf')
                min_2 = float('inf')
                arg_min_1 = float('inf')
                arg_min_2 = float('inf')
                for x_j in N:
                    if (x_i!=x_j):
                        value = g[x_i][x_j];
                        if (value<min_1):
                            min_2 = min_1
                            min_1 = value
                            arg_min_2 = arg_min_1
                            arg_min_1 = x_j
                        elif (value<min_2):
                            min_2 = value
                            arg_min_2 = x_j

                p[l][x_i] = arg_min_1
                f[l][x_i] = g[x_i][arg_min_1]
                phi[l][x_i] = g[x_i][arg_min_2]
                coord = g_type[x_i][arg_min_1]
                coord_2 = g_type[x_i][arg_min_2]
                q_route[l][x_i] = (q_route[coord[1]][arg_min_1] + [x_i]) if (coord[0] == 0) else (q_route_2[coord[1]][arg_min_1] + [x_i])
                q_route_2[l][x_i] = (q_route[coord_2[1]][arg_min_2] + [x_i]) if (coord_2[0] == 0) else (q_route_2[coord_2[1]][arg_min_2] + [x_i])
    return (f, phi, p, q_route, q_route_2)



# Combines paths to construct routes
def construct_q_routes(h,truck_capacity,N,distance,values,values_pos,quantities):
    f_l, phi_l, p_l, q_route_l, q_route_2_l = construct_q_paths(h,truck_capacity,N,distance,values,values_pos,quantities,'left')
    f_r, phi_r, p_r, q_route_r, q_route_2_r = construct_q_paths(h,truck_capacity,N,distance,values,values_pos,quantities,'right')

    # Initialize the routes
    psi = {}; psi_route = {};
    for l in range(len(values)):
        psi[l] = {}; psi_route[l] = {};
        for n in N:
            psi[l][n] = {}; psi_route[l][n] = {};

    for l in range(len(values)):
        for n in N:
            min_w = quantities[n]
            max_w = (values[l] + quantities[n])
            min_route = []
            min_val = float('inf')
            for l_1 in range(len(values)):
                q = values[l_1]
                if q>=min_w and q<max_w:
                    l_2 = values_pos[values[l] + quantities[n]-q]
                    if p_l[l_1][n]!=p_r[l_2][n] or (p_l[l_1][n] == p_r[l_2][n] and p_l[l_1][n]==h):
                        val = f_l[l_1][n]+f_r[l_2][n]
                        route = q_route_l[l_1][n]+list(reversed(q_route_r[l_2][n]))[1:]
                    else:
                        if f_l[l_1][n]+phi_r[l_2][n]<phi_l[l_1][n]+f_r[l_2][n]:
                            val = f_l[l_1][n]+phi_r[l_2][n]
                            route = q_route_l[l_1][n]+list(reversed(q_route_2_r[l_2][n]))[1:]
                        else:
                            val = phi_l[l_1][n]+f_r[l_2][n]
                            route = q_route_2_l[l_1][n]+list(reversed(q_route_r[l_2][n]))[1:]
                    if val<min_val:
                        min_val = val
                        min_route = route
            psi[l][n] = min_val
            psi_route[l][n] = min_route
    return (psi, psi_route)

# Combines paths to construct routes
def construct_q_routes_debug(h,truck_capacity,N,distance,values,values_pos,quantities):
    f_l, phi_l, p_l, q_route_l, q_route_2_l = construct_q_paths_debug(h,truck_capacity,N,distance,values,values_pos,quantities,'left')
    f_r, phi_r, p_r, q_route_r, q_route_2_r = construct_q_paths_debug(h,truck_capacity,N,distance,values,values_pos,quantities,'right')

    # Initialize the routes
    psi = {}; psi_route = {};
    for l in range(len(values)):
        psi[l] = {}; psi_route[l] = {};
        for n in N:
            psi[l][n] = {}; psi_route[l][n] = {};

    for l in range(len(values)):
        for n in N:
            min_w = quantities[n]
            max_w = (values[l] + quantities[n])
            min_route = []
            min_val = float('inf')
            min_coord = ()
            for l_1 in range(len(values)):
                q = values[l_1]
                if q>=min_w and q<max_w:
                    l_2 = values_pos[values[l] + quantities[n]-q]
                    if p_l[l_1][n]!=p_r[l_2][n] or (p_l[l_1][n] == p_r[l_2][n] and p_l[l_1][n]==h):
                        val = f_l[l_1][n]+f_r[l_2][n]
                        coord = (0,l_1,l_2)
                    else:
                        if f_l[l_1][n]+phi_r[l_2][n]<phi_l[l_1][n]+f_r[l_2][n]:
                            val = f_l[l_1][n]+phi_r[l_2][n]
                            coord = (1,l_1,l_2)
                        else:
                            val = phi_l[l_1][n]+f_r[l_2][n]
                            coord = (2,l_1,l_2)
                    if val<min_val:
                        min_val = val
                        min_coord = coord
            psi[l][n] = min_val
            if len(min_coord)>0:
                if min_coord[0] == 0:
                    min_route = q_route_l[min_coord[1]][n]+list(reversed(q_route_r[min_coord[2]][n]))[1:]
                elif min_coord[0] == 1:
                    min_route = q_route_l[min_coord[1]][n]+list(reversed(q_route_2_r[min_coord[2]][n]))[1:]
                elif min_coord[0] == 2:
                    min_route = q_route_2_l[min_coord[1]][n]+list(reversed(q_route_r[min_coord[2]][n]))[1:]
            psi_route[l][n] = min_route
    return (psi, psi_route)



# Calculate lower bound of the optimization problem, along with auxiliary quantities
def lower_bound(H,capacities,N,quantities,distance,mu,lamb):
    # Calculate the bound of the problem and the routes that help us reach that bound
    psi = {};
    psi_route = {};
    b = {};
    W = {}
    for h in H:
        truck_capacity = capacities[h]
        values = possible_values(quantities,capacities[h])
        W[h] = values
        values_pos = dict(zip(values,list(range(len(values)))))
        psi_h, psi_route_h = construct_q_routes_debug(h,truck_capacity,N,distance,values,values_pos,quantities)
        psi[h] = psi_h
        psi_route[h] = psi_route_h
        b[h] = {}
        for l in range(len(values)):
            b[h][l] = {}
            for n in N:
                b[h][l][n] = (psi[h][l][n]-mu[h])*quantities[n]/values[l]

    min_b = {}
    arg_l = {}
    for h in H:
        min_b[h] = {}
        arg_l[h] = {}
        for n in N:
            min_b[h][n] = np.min(np.array([b[h][l][n] for l in range(len(W[h]))]))
            arg_l[h][n] = np.argmin(np.array([b[h][l][n] for l in range(len(W[h]))]))

    min_min_b = {}
    arg_h = {}
    arg_arg_l = {}
    arg_route = {}
    for n in N:
        i = True
        min_min_b[n] = np.min(np.array([min_b[h][n] for h in H]))
        arg_h[n] = H[np.argmin(np.array([min_b[h][n] for h in H]))]
        arg_arg_l[n] = arg_l[arg_h[n]][n]
        arg_route[n] = psi_route[arg_h[n]][arg_l[arg_h[n]][n]][n]

    # Compute theta and rho
    # First we calculate the number of times each node is visited
    visits = {}
    for n in N:
        visits[n] = {}
        for n_2 in N:
            visits[n][n_2] = np.sum([n_2 == i for i in arg_route[n]])

    theta = {}
    for j in N:
        theta[j] = 0
        for i in N:
            theta[j] = theta[j] + (quantities[i]+0.0)/W[arg_h[i]][arg_arg_l[i]]*visits[i][j]

    rho = {}
    for h in H:
        rho[h] = 0
    for i in N:
        rho[arg_h[i]] = rho[arg_h[i]] + (quantities[i]+0.0)/W[arg_h[i]][arg_arg_l[i]]

    # Construct the duality vectors
    u = {}
    for n in N:
        u[n] = min_min_b[n] + lamb[n]

    z_lb = np.sum([u[n] for n in N]) + np.sum([mu[h] for h in H])

    return (z_lb, theta, rho, u)



def optimize_lower_bound(iterations, z_ub, epsilon, H,capacities,N,quantities,distance_dictionary):

    # Initialize the parameters
    lamb = {}
    mu = {}
    for j in N:
        lamb[j] = 0
    for h in H:
        mu[h] = 0
    max_val = -float('inf')
    values = []

    for i in range(iterations):
        distance = reduced_cost_dict(lamb, distance_dictionary, N)

        z_lb, theta, rho, u = lower_bound(H,capacities,N,quantities,distance,mu,lamb)
        values.append(z_lb)
        if z_lb > max_val:
            max_val = copy.deepcopy(z_lb)
            u_opt = copy.deepcopy(u)
            v_opt = copy.deepcopy(mu)
            lamb_opt = copy.deepcopy(lamb)

        # Compute the new parameters
        gamma = (z_ub - z_lb)/(np.sum((np.array(theta.values())-1)**2) + np.sum((np.array(rho.values())-1)**2))
        # New lambda
        for j in N:
            lamb[j] = lamb[j] - epsilon*gamma*(theta[j]-1)
        for h in H:
            mu[h] = np.min([mu[h] - epsilon*gamma*(rho[h]-1),0])
        print(z_lb)

        if np.sum(np.abs(lamb.values())) > 10**16:
            raise ValueError('Lambda exploding')

        # Rule for updating epsilon
        if len(values)>=7:
            grad = [values[i+1]-values[i] for i in range(len(values)-1)]
            jumps = [np.sign(grad[i+1])!=np.sign(grad[i]) for i in range(len(grad)-1)]
            if np.sum(jumps[len(jumps)-5:len(jumps)])>=3:
                epsilon = epsilon/1.5
                print ('new epsilon:%f' % epsilon)
                values = []
            if np.sum(np.array(grad[len(grad)-5:len(grad)])>0) >= 5:

                epsilon = epsilon*1.2
                print ('new epsilon:%f' % epsilon)
                values = []
        if (np.array_equal(np.array(theta.values()), np.ones(len(N))) and np.array_equal(np.array(rho.values()), np.ones(len(H)))):
            print("reached zero gradient")
            return (max_val,u_opt,v_opt,lamb_opt)
    return (max_val,u_opt,v_opt,lamb_opt)

def optimize_lower_bound_c(iterations, z_ub, epsilon, H,capacities,N,quantities,distance_mat):

    # Initialize the parameters
    mu_ = np.zeros(len(H), dtype = "float64")
    lamb_ = np.zeros(len(N), dtype = "float64")
    max_val = -float('inf')
    values = []

    H_ = (np.array(range(len(H)))+len(N)).astype(int)
    N_ = np.array(range(len(N))).astype(int)
    capacities_ = np.array([capacities[h] for h in H]).astype(int)
    quantities_ = np.array([quantities[n] for n in N]).astype(int)

    for i in range(iterations):
        distance_ = reduced_cost_mat(lamb_, distance_mat, N_)
        results_c = cpp_lower_bounds.lower_bound(H_,capacities_,N_,quantities_,distance_,mu_,lamb_)
        z_lb = results_c["z_lb"]
        theta = np.array(results_c["theta"]).astype("float64")
        rho = np.array(results_c["rho"]).astype("float64")
        u = np.array(results_c["u"]).astype("float64")
        print(z_lb)
        values.append(z_lb)
        if z_lb > max_val:
            max_val = copy.deepcopy(z_lb)
            u_opt = copy.deepcopy(u)
            v_opt = copy.deepcopy(mu_)
            lamb_opt = copy.deepcopy(lamb_)

        # Compute the new parameters
        gamma = (z_ub - z_lb)/(np.sum((theta-1)**2) + np.sum((rho-1)**2))
        # New lambda
        lamb_ = np.array(lamb_ - epsilon*gamma*(theta-1))
        # New mu
        mu_ = np.array(np.minimum(mu_ - epsilon*gamma*(rho-1),0)).astype('float64')

        if np.sum(np.abs(lamb_)) > 10**16:
            raise ValueError('Lambda exploding')

        # Rule for updating epsilon
        if len(values)>=7:
            grad = [values[i+1]-values[i] for i in range(len(values)-1)]
            jumps = [np.sign(grad[i+1])!=np.sign(grad[i]) for i in range(len(grad)-1)]
            if np.sum(jumps[len(jumps)-5:len(jumps)])>=3:
                epsilon = epsilon/1.5
                print ('new epsilon:%f' % epsilon)
                values = []
            if np.sum(np.array(grad[len(grad)-5:len(grad)])>0) >= 5:

                epsilon = epsilon*1.2
                print ('new epsilon:%f' % epsilon)
                values = []
        if (np.array_equal(theta, np.ones(len(N_))) and np.array_equal(rho, np.ones(len(H_)))):
            print("reached zero gradient")
            return (max_val,u_opt,v_opt,lamb_opt)

    return (max_val,u_opt,v_opt,lamb_opt)
