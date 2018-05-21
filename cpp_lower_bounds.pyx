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
   QPaths construct_q_paths_(
      int h,
      int truck_capacity,
      vector[int] N,
      vector[vector[double]] distance,
      vector[int] values,
      map[double,int] values_pos,
      vector[int] quantities,
      string direction
      )
   QRoutes construct_q_routes_(
      int h,
      int truck_capacity,
      vector[int] N,
      vector[vector[double]] distance,
      vector[int] values,
      map[double,int] values_pos,
      vector[int] quantities
      )

cpdef construct_q_paths(h_,truck_capacity_,N_,distance_,values_,values_pos_,quantities_,direction_):
    cdef:
        int h = h_
        int truck_capacity = truck_capacity_
        vector[int] N = N_
        vector[vector[double]] distance = distance_
        vector[int] values = values_
        map[double,int] values_pos = values_pos_
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
      map[double,int] values_pos = values_pos_
      vector[int] quantities = quantities_

   cdef QRoutes qroutes = construct_q_routes_(h,truck_capacity,N,distance,values,values_pos,quantities)
   return qroutes
