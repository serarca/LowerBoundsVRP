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


using namespace std;

int main(){

   int h = 3;
   int n = 10;
   vector<int> capacities(h, 5);
   vector<int> quantities(n, 1);
   vector<int> N(n,0);
   for (int i = 0; i<n; i++){
      N[i] = i;
   }
   vector<int> H(h,0);
   for (int i = 0; i<h; i++){
      H[i] = n+i;
   }
   vector<double> mu(h,0);
   vector<double> lamb(n,0);

   vector<double> x(n+h);
   vector<double> y(n+h);


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


   vector<vector<double>> distance_dict(n+h, vector<double>(n+h));
   for (int i = 0; i < n+h; i++){
      for (int j = 0; j < n+h; j++){
         distance_dict[i][j] = sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2));
      }
   }



   LowerBound lb = lower_bound_(H, capacities, N, quantities, distance_dict, mu, lamb);
   cout<<lb.z_lb<<endl;


   return 1;
}
