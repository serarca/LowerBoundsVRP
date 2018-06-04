#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <limits>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <random>
#include <functional>
#include<cmath>
#include "lower_bounds.h"
#include "baldacci.h"
#include "lower_bounds.cpp"
#include "baldacci.cpp"
template class std::set<tuple<int,int>>;
# include <chrono>
using  ms = chrono::milliseconds;
using get_time = chrono::steady_clock ;
#include "prettyprint.hpp"


void print_set(set<tuple<int,int>> s){
   for (auto n:s)
      cout<<"("<<get<0>(n)<<","<<get<1>(n)<<")";
   cout<<endl;
}

void print_path(Path p){
   cout<<"\t Path: ";
   for (auto n:p.path){
      cout<<n<<" ";
   }
   cout<<"Set: ";
   for (auto n:p.nodes){
      cout<<n<<" ";
   }
   cout<<"Cost: "<<p.cost<<" ";
   cout<<"Load: "<<p.load<<" ";
   cout<<"Lower Bound: "<<p.lower_bound<<" "<<endl;
}

void print_route(Route r){
   cout<<"\t Route_l: ";
   for (auto n:(r.path_l->path)){
      cout<<n<<" ";
   }
   cout<<"\t Route_r: ";
   for (auto n:(r.path_r->path)){
      cout<<n<<" ";
   }
   cout<<"Set: ";
   for (auto n:r.nodes){
      cout<<n<<" ";
   }
   cout<<"Cost: "<<r.cost<<" ";
   cout<<"Load: "<<r.load<<" ";
   cout<<"Median: "<<r.median<<" "<<endl;
}

void print_sroute(SimpleRoute r){
   cout<<"\t Route: ";
   for (auto n:(r.path)){
      cout<<n<<" ";
   }
   cout<<"Cost: "<<r.cost<<" ";
   cout<<"Geo Cost: "<<r.geo_cost<<" ";
   cout<<"Load: "<<r.load<<" ";
   cout<<"Median: "<<r.median<<" "<<endl;
}

void print_Paths(vector<list<Path>> paths){
   int i = 0;
   for (auto end:paths){
      cout<<"End Node: "<<i;
      cout<<"Length: "<<end.size()<<endl;
      i+=1;

      for (auto p:end){
         print_path(p);
      }
   }
}

void print_Routes(vector<list<Route>> routes){
   int i = 0;
   for (auto end:routes){
      cout<<"End Node: "<<i;
      cout<<"Length: "<<end.size()<<endl;
      i+=1;

      for (auto r:end){
         print_route(r);
      }
   }
}

void print_sRoutes(list<SimpleRoute> end){
      cout<<"Length: "<<end.size()<<endl;

      for (auto r:end){
         print_sroute(r);
      }
}

void print_SRoutes(vector<list<SimpleRoute>> routes){
   int i = 0;
   for (auto end:routes){
      cout<<"End Node: "<<i;
      cout<<"Length: "<<end.size()<<endl;
      i+=1;

      for (auto r:end){
         print_sroute(r);
      }
   }
}

void p_v(vector<double> vec){
   cout<<vec<<endl;
}

void p_v_v(vector<vector<double>> vec){
   for (int i = 0; i< (int)vec.size(); i++){
      cout<<i<<":"<<vec[i]<<endl;
   }

}




int main(){

   int H_len = 3;
   int N_len = 3;
   vector<int> capacities(H_len, 10);
   vector<int> quantities(N_len, 1);
   vector<int> N(N_len,0);
   for (int i = 0; i<N_len; i++){
      N[i] = i;
   }
   vector<int> H(H_len,0);
   for (int i = 0; i<H_len; i++){
      H[i] = N_len+i;
   }
   vector<double> mu(H_len,0);
   vector<double> lamb(N_len,0);

   vector<double> x(N_len+H_len);
   vector<double> y(N_len+H_len);


   // First create an instance of an engine.
   random_device rnd_device;
   // Specify the engine and distribution.
   mt19937 mersenne_engine(rnd_device());
   mersenne_engine.seed(std::random_device()());
   uniform_real_distribution<double> dist(0, 1);

   auto gen = std::bind(dist, mersenne_engine);
   generate(begin(x), end(x), gen);

   mersenne_engine.seed(std::random_device()());
   auto gen2 = std::bind(dist, mersenne_engine);
   generate(begin(y), end(y), gen2);


   vector<vector<double>> geo_distance(N_len+H_len, vector<double>(N_len+H_len));
   for (int i = 0; i < N_len+H_len; i++){
      for (int j = 0; j < N_len+H_len; j++){
         geo_distance[i][j] = sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2));
      }
   }

   //LowerBound lb = lower_bound_(H, capacities, N, quantities, distance_dict, mu, lamb);
   //cout<<lb.z_lb<<endl;

   /*
   int Delta = 700;
   double gamma = 12;
   int h = H[0];
   int capacity = capacities[0];
   string direction = "left";
   //vector<list<Path>> p = GENPATH(Delta, gamma, h, capacity, N, quantities, distance_dict, direction);
   auto start = get_time::now();
   vector<list<SimpleRoute>> r = GENROUTE(Delta, gamma, h, capacity, N, quantities, distance_dict);
   auto end = get_time::now();
   auto diff = end - start;
   cout<<"Elapsed time is :  "<< chrono::duration_cast<ms>(diff).count()<<" ms "<<endl;
   print_SRoutes(r);
   */

   // Debugging first lower bounds
   double z_ub = 20;
   int iterations = 100;
   double epsilon = 0.1;
   DualSolution sol = lower_bound_optimizer_M1(iterations, z_ub, epsilon, H, capacities, N, quantities,geo_distance);



   // Debugging the lower_bound_M2
   int Delta = 1;
   int Delta_zero = 10;
   double gamma = 5;
   double gamma_zero = 5;
   int len_N = N.size();
   int len_H = H.size();
   int sub_iterations = 100;

   // Define the vector of routes for each truck
   vector<list<SimpleRoute>> Routes(len_H);

   // Calculate reduced costs
   vector<vector<double>> distance_dict = reduced_cost_matrix(geo_distance, sol.lamb, sol.v);

   // We start by generating routes for all of the trucks
   for (auto h:H){
      Routes[h - len_N] = GENROUTE(Delta, gamma, h, capacities[h - len_N], N, quantities, distance_dict, geo_distance);
   }

   while(true){
      // We pass these routes to the algorithm that calculates the lower bound
      DualSolution ds = lower_bound_optimizer_M2(sub_iterations, Routes, z_ub, epsilon, H, capacities, N, quantities, sol.v, sol.lamb);


      vector<list<SimpleRoute>> newRoutes(len_H);
      vector<vector<double>> new_distance_dict = reduced_cost_matrix(geo_distance, ds.lamb, ds.v);
      int new_routes_count = 0;
      for (auto h:H){
         newRoutes[h - len_N] = GENROUTE(Delta_zero, gamma_zero, h, capacities[h - len_N], N, quantities, new_distance_dict, geo_distance);
         for (auto route:newRoutes[h - len_N]){
            if (route.cost < 0){
               Routes[h - len_N].push_back(route);
               ++new_routes_count;
            }
         }
      }
      cout<<"Negative routes: "<<new_routes_count<<endl;
      if (new_routes_count == 0){
         break;
      }

   }

   return 1;
}
