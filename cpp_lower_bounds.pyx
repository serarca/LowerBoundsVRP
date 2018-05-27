import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.string cimport string

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
