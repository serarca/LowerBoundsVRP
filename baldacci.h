#pragma once

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <list>
#include <set>

using namespace std;

//The struct of paths
struct Path{
   list<int> path;
   set<int> nodes;
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
   set<int> nodes;
   double cost;
   int load;
   int median;
};

//The struct of paths
struct SimpleRoute{
   list<int> path;
   int index_l;
   int index_r;
   double cost;
   int load;
   int median;
   double geo_cost;
   int truck;
};



vector<list<Path>> GENPATH(
   int Delta,
   double gamma,
   int h,
   int capacity,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> distance_dict,
   string direction);

list<SimpleRoute> GENROUTE(
   int Delta,
   double gamma,
   int h,
   int capacity,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> distance_dict,
   vector<vector<double>> geo_distance);

void optimize_lower_bound(
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
