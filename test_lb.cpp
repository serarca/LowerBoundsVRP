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




int main(){

   int H_len = 2;
   int N_len = 10;
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

   int h = H[0];
   int truck_capacity = capacities[0];
   vector<vector<double>> distance = geo_distance;
   distance[0][1] = -distance[0][1]*100;

   PossibleValues p = possible_values(quantities, truck_capacity);
   string direction = "left";

   LB_GENPATH lb = path_lower_bound(
      h,
      truck_capacity,
      N,
      distance,
      p.values,
      p.values_pos,
      quantities,
      direction
   );

   p_v_v(lb.F);

   return 1;
}
