#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <limits>
#include <iterator>



using namespace std;

// The struct of results
struct QPaths {
   vector<vector<double>> f;
   vector<vector<double>> phi;
   vector<vector<int>> p;
   vector<vector<vector<int>>> q_route;
   vector<vector<vector<int>>> q_route_2;
};


QPaths construct_q_paths_(
   int h,
   int truck_capacity,
   vector<int> N,
   vector<vector<double>> distance_dict,
   vector<int> values,
   map<double,int> values_pos,
   vector<int> quantities,
   string direction
){
   //Calculate lengths
   int len_values = values.size();
   int len_N = N.size();

   //Construct infinity
   double inf = numeric_limits<double>::infinity();

   //Initialize the routes
   vector<vector<double>> f(len_values, vector<double> (len_N));
   vector<vector<double>> phi(len_values, vector<double> (len_N));
   vector<vector<int>> p(len_values, vector<int> (len_N));
   vector<vector<vector<int>>> q_route(len_values, vector<vector<int>> (len_N, vector<int>(0)));
   vector<vector<vector<int>>> q_route_2(len_values, vector<vector<int>> (len_N, vector<int>(0)));

   for(int l = 0; l < len_values; l++) {
      for (int n = 0; n < len_N; n++) {
         f[l][n] = inf;
         phi[l][n] = inf;
      }
   }
   // Initialize the routes
   for (int n = 0; n < len_N; n++) {
      int q = quantities[n];
      if (q <= truck_capacity){
         int l = values_pos[q];
         if (direction == "left"){
            f[l][n] = distance_dict[h][N[n]];
         } else {
            f[l][n]  = distance_dict[N[n]][h];
         }
         p[l][n] = h;
         const int args[] = {h, N[n]};
         q_route[l][n].insert(q_route[l][n].end(), args, args+2);
      }
   }


   for(int l = 0; l < len_values; l++) {
      int Q = values[l];
      vector<vector<double>> g(len_N,vector<double>(len_N));
      vector<vector<vector<int>>> g_type(len_N,vector<vector<int>>(len_N,vector<int>(2,0)));
      for (int x_i = 0; x_i < len_N; x_i++){
         int q_p = Q - quantities[x_i];
         if (q_p > 0){
            int l_p = values_pos[q_p];
            for (int x_j = 0; x_j < len_N; x_j++){
               if (x_i != x_j){
                  if (p[l_p][x_j]!=x_i){
                     g[x_i][x_j] = f[l_p][x_j] + distance_dict[N[x_j]][N[x_i]];
                     // We save a boolean and the quantity at which g is calc
                     g_type[x_i][x_j][0] = 0;
                     g_type[x_i][x_j][1] = l_p;
                  } else {
                     g[x_i][x_j] = phi[l_p][x_j] + distance_dict[N[x_j]][N[x_i]];
                     g_type[x_i][x_j][0] = 1;
                     g_type[x_i][x_j][1] = l_p;
                  }
               }
            }
         }
      }
      for (int x_i = 0; x_i < len_N; x_i++){
         int q_p = Q - quantities[x_i];
         if (q_p > 0){
            // Find the minimum
            int arg_min_1;
            int arg_min_2;
            double min_1 = inf;
            double min_2 = inf;
            for (int x_j = 0; x_j < len_N; x_j++){
               if (x_i!=x_j){
                  double value = g[x_i][x_j];
                  if (value<min_1){
                     min_2 = min_1;
                     min_1 = value;
                     arg_min_2 = arg_min_1;
                     arg_min_1 = x_j;
                  } else if (value<min_2) {
                     min_2 = value;
                     arg_min_2 = x_j;
                  }
               }
            }
            p[l][x_i] = arg_min_1;
            f[l][x_i] = min_1;
            phi[l][x_i] = min_2;
            vector<int> &coord = g_type[x_i][arg_min_1];
            vector<int> &coord_2 = g_type[x_i][arg_min_2];

            q_route[l][x_i] = (coord[0] == 0) ? q_route[coord[1]][arg_min_1] : q_route_2[coord[1]][arg_min_1];
            q_route_2[l][x_i] = (coord_2[0] == 0) ? q_route[coord_2[1]][arg_min_2] : q_route_2[coord_2[1]][arg_min_2];
            q_route[l][x_i].push_back(x_i);
            q_route_2[l][x_i].push_back(x_i);

         }
      }
   }


   QPaths qpaths;
   qpaths.f = f;
   qpaths.phi = phi;
   qpaths.p = p;
   qpaths.q_route = q_route;
   qpaths.q_route_2 = q_route_2;

   return qpaths;

}
