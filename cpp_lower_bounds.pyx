import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.string cimport string

cdef extern from "lower_bounds.h":
   struct QPaths:
      map[int,map[string,double]] f
      map[int,map[string,double]] phi
      map[int,map[string,string]] p
      map[int,map[string,vector[string]]] q_route
      map[int,map[string,vector[string]]] q_route_2
   QPaths construct_q_paths_(
      string h,
      int truck_capacity,
      vector[string] N,
      map[string,map[string,double]] distance,
      vector[int] values,
      map[double,int] values_pos,
      map[string,int] quantities,
      string direction
      );

cpdef construct_q_paths(h_,truck_capacity_,N_,distance_,values_,values_pos_,quantities_,direction_):
    cdef:
        string h = h_
        int truck_capacity = truck_capacity_
        vector[string] N = N_
        map[string,map[string,double]] distance = distance_
        vector[int] values = values_
        map[double,int] values_pos = values_pos_
        map[string,int] quantities = quantities_
        string direction = direction_

    cdef QPaths qpaths = construct_q_paths_(h,truck_capacity,N,distance,values,values_pos,quantities,direction)
    return qpaths
