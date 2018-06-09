#pragma once

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <list>
#include <set>
#include "lower_bounds.h"
#include <bitset>

using namespace std;

typedef std::bitset<100> bit_set;


//The struct of paths
struct Path{
   list<int> path;
   bit_set nodes;
   double cost;
   double lower_bound;
   int load;
   int end;
};

//The struct of paths
struct Route{
   list<Path>::iterator path_l;
   list<Path>::iterator path_r;
   int index_l;
   int index_r;
   bit_set nodes;
   double cost;
   int load;
   int median;
};




vector<list<Path>> GENPATH(
   int Delta,
   double gamma,
   int h,
   int capacity,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> distance_dict,
   string direction,
   bool &terminated
);

list<SimpleRoute> GENROUTE(
   int Delta,
   double gamma,
   int h,
   int capacity,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> distance_dict,
   vector<vector<double>> geo_distance,
   bool &terminated
);

DualSolution optimize_lower_bound_M2(
   int sub_iterations,
   double z_ub,
   int Delta,
   int Delta_zero,
   double gamma,
   double gamma_zero,
   double epsilon,
   vector<int> H,
   vector<int> capacities,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> geo_distance,
   vector<double> mu,
   vector<double> lamb
);

DualSolution construct_lower_bound(
   int iterations_grad_m1,
   int iterations_grad_m2,
   int iterations_m2,
   double z_ub,
   int Delta,
   int Delta_zero,
   int Delta_final,
   double gamma,
   double gamma_zero,
   double gamma_final,
   double epsilon,
   vector<int> H,
   vector<int> capacities,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> geo_distance
);
