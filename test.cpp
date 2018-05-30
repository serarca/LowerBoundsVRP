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
#include "lower_bounds.cpp"
#include "baldacci.cpp"
template class std::set<tuple<int,int>>;
# include <chrono>
using  ms = chrono::milliseconds;
using get_time = chrono::steady_clock ;



using namespace std;
void print_set(set<tuple<int,int>> s){
   for (auto n:s)
      cout<<"("<<get<0>(n)<<","<<get<1>(n)<<")";
   cout<<endl;
}


int main(){

   int H_len = 30;
   int N_len = 200;
   vector<int> capacities(H_len, 5);
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


   vector<vector<double>> distance_dict(N_len+H_len, vector<double>(N_len+H_len));
   for (int i = 0; i < N_len+H_len; i++){
      for (int j = 0; j < N_len+H_len; j++){
         distance_dict[i][j] = sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2));
      }
   }

   auto start = get_time::now();
   LowerBound lb = lower_bound_(H, capacities, N, quantities, distance_dict, mu, lamb);
   cout<<lb.z_lb<<endl;
   auto end = get_time::now();
   auto diff = end - start;
   cout<<"Elapsed time is :  "<< chrono::duration_cast<ms>(diff).count()<<" ms "<<endl;
   
   /*
   int Delta = 200;
   double gamma = 12;
   int h = H[0];
   int capacity = capacities[0];
   string direction = "left";
   vector<list<Path>> p = GENPATH(Delta, gamma, h, capacity, N, quantities, distance_dict, direction);
   auto start = get_time::now();
   vector<list<Route>> r = GENROUTE(Delta, gamma, h, capacity, N, quantities, distance_dict);
   auto end = get_time::now();
   auto diff = end - start;
   cout<<"Elapsed time is :  "<< chrono::duration_cast<ms>(diff).count()<<" ms "<<endl;
   */
   return 1;
}
