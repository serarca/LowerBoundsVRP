import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.list cimport list
from libcpp.set cimport set



cdef extern from "lower_bounds.h":
   struct QPaths:
      vector[vector[double]] f
      vector[vector[double]] phi
      vector[vector[int]] p
      vector[vector[vector[int]]] q_route
      vector[vector[vector[int]]] q_route_2
   struct QRoutes:
      vector[vector[double]] psi
      vector[vector[vector[int]]] psi_route
   struct LowerBound:
      double z_lb
      vector[double] theta
      vector[double] rho
      vector[double] u
   struct DualSolution:
      double z_lb
      vector[double] lamb
      vector[double] u
      vector[double] v

   QPaths construct_q_paths_(
      int h,
      int truck_capacity,
      vector[int] N,
      vector[vector[double]] distance,
      vector[int] values,
      map[int,int] values_pos,
      vector[int] quantities,
      string direction
      )
   QRoutes construct_q_routes_(
      int h,
      int truck_capacity,
      vector[int] N,
      vector[vector[double]] distance,
      vector[int] values,
      map[int,int] values_pos,
      vector[int] quantities
      )

   LowerBound lower_bound_(
      vector[int] H,
      vector[int] capacities,
      vector[int] N,
      vector[int] quantities,
      vector[vector[double]] distance_dict,
      vector[double] mu,
      vector[double] lamb
      )
   DualSolution lower_bound_optimizer_M1(
      int iterations,
      double z_ub,
      double epsilon,
      vector[int] H,
      vector[int] capacities,
      vector[int] N,
      vector[int] quantities,
      vector[vector[double]] geo_distance);


cdef extern from "baldacci.h":
   struct Path:
      list[int] path
      set[int] nodes
      double cost
      double lower_bound
      int load
      int end

   struct SimpleRoute:
      list[int] path
      int index_l
      int index_r
      double cost
      int load
      int median
      double geo_cost
      int truck

   vector[list[Path]] GENPATH(
      int Delta,
      double gamma,
      int h,
      int capacity,
      vector[int] N,
      vector[int] quantities,
      vector[vector[double]] distance_dict,
      string direction
      )

   list[SimpleRoute] GENROUTE(
      int Delta,
      double gamma,
      int h,
      int capacity,
      vector[int] N,
      vector[int] quantities,
      vector[vector[double]] distance_dict,
      vector[vector[double]] geo_distance
      )

cpdef construct_q_paths(h_,truck_capacity_,N_,distance_,values_,values_pos_,quantities_,direction_):
    cdef:
        int h = h_
        int truck_capacity = truck_capacity_
        vector[int] N = N_
        vector[vector[double]] distance = distance_
        vector[int] values = values_
        map[int,int] values_pos = values_pos_
        vector[int] quantities = quantities_
        string direction = direction_

    cdef QPaths qpaths = construct_q_paths_(h,truck_capacity,N,distance,values,values_pos,quantities,direction)
    return qpaths

cpdef construct_q_routes(h_,truck_capacity_,N_,distance_,values_,values_pos_,quantities_):
   cdef:
      int h = h_
      int truck_capacity = truck_capacity_
      vector[int] N = N_
      vector[vector[double]] distance = distance_
      vector[int] values = values_
      map[int,int] values_pos = values_pos_
      vector[int] quantities = quantities_

   cdef QRoutes qroutes = construct_q_routes_(h,truck_capacity,N,distance,values,values_pos,quantities)
   return qroutes

cpdef lower_bound(H_, capacities_, N_, quantities_, distance_, mu_, lamb_):
   cdef:
      vector[int] H = H_,
      vector[int] capacities = capacities_,
      vector[int] N = N_,
      vector[int] quantities = quantities_,
      vector[vector[double]] distance = distance_,
      vector[double] mu = mu_,
      vector[double] lamb = lamb_
   cdef LowerBound lb = lower_bound_(H, capacities, N, quantities, distance, mu, lamb)
   return lb

cpdef GENPATH_(Delta_, gamma_, h_, capacity_, N_, quantities_, distance_dict_, direction_):
   cdef:
      int Delta = Delta_
      double gamma = gamma_
      int h = h_
      int capacity = capacity_
      vector[int] N = N_
      vector[int] quantities = quantities_
      vector[vector[double]] distance_dict = distance_dict_
      string direction = direction_
   cdef vector[list[Path]] paths = GENPATH(Delta, gamma, h, capacity, N, quantities, distance_dict, direction)
   return paths

cpdef GENROUTE_(Delta_, gamma_, h_, capacity_, N_, quantities_, distance_dict_, geo_distance_):
   cdef:
      int Delta = Delta_
      double gamma = gamma_
      int h = h_
      int capacity = capacity_
      vector[int] N = N_
      vector[int] quantities = quantities_
      vector[vector[double]] distance_dict = distance_dict_
      vector[vector[double]] geo_distance = geo_distance_
   cdef list[SimpleRoute] routes = GENROUTE(Delta, gamma, h, capacity, N, quantities, distance_dict, geo_distance)
   return routes

cpdef lower_bound_optimizer_M1_(iterations_, z_ub_, epsilon_, H_, capacities_, N_, quantities_, geo_distance_):
   cdef:
      int iterations = iterations_
      double z_ub = z_ub_
      double epsilon = epsilon_
      vector[int] H = H_
      vector[int] capacities = capacities_
      vector[int] N = N_
      vector[int] quantities = quantities_
      vector[vector[double]] geo_distance = geo_distance_
   cdef DualSolution dual_solution = lower_bound_optimizer_M1(iterations, z_ub, epsilon, H, capacities, N, quantities, geo_distance)
   return dual_solution
