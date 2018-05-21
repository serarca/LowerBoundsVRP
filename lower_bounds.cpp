#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <limits>
#include <iterator>



using namespace std;

// The struct of results
struct QPaths {
   map<int,map<string,double>> f;
   map<int,map<string,double>> phi;
   map<int,map<string,string>> p;
   map<int,map<string,vector<string>>> q_route;
   map<int,map<string,vector<string>>> q_route_2;
};


QPaths construct_q_paths_(
   string h,
   int truck_capacity,
   vector<string> N,
   map<string,map<string,double>> distance_dict,
   vector<int> values,
   map<double,int> values_pos,
   map<string,int> quantities,
   string direction
){
   //Calculate lengths
   int len_values = values.size();
   int len_N = N.size();

   //Construct infinity
   double inf = numeric_limits<double>::infinity();

   //Initialize the routes
   map<int,map<string,double>> f;
   map<int,map<string,double>> phi;
   map<int,map<string,string>> p;
   map<int,map<string,vector<string>>> q_route;
   map<int,map<string,vector<string>>> q_route_2;
   for(int l = 0; l < len_values; l++) {
      for (vector<string>::iterator n_it = N.begin(); n_it != N.end(); ++n_it) {
         string n = *n_it;
         f[l][n] = inf;
         phi[l][n] = inf;
         p[l][n] = inf;
      }
   }

   // Initialize the routes
   for (vector<string>::iterator n_it = N.begin(); n_it != N.end(); ++n_it) {
      string n = *n_it;
      int q = quantities[n];
      if (q <= truck_capacity){
         int l = values_pos[q];
         if (direction == "left"){
            f[l][n] = distance_dict[h][n];
         } else {
            f[l][n]  = distance_dict[n][h];
         }
         p[l][n] = h;
         const string args[] = {h, n};
         q_route[l][n].insert(q_route[l][n].end(), args, args+2);
      }
   }


   for(int l = 0; l < len_values; l++) {
      int Q = values[l];
      map<string,map<string,double>> g;
      map<string,map<string,vector<string>>> g_type;
      for (vector<string>::iterator x_i_it = N.begin(); x_i_it != N.end(); ++x_i_it){
         string x_i = *x_i_it;
         int q_p = Q - quantities[x_i];
         if (q_p > 0){
            int l_p = values_pos[q_p];
            for (vector<string>::iterator x_j_it = N.begin(); x_j_it != N.end(); ++x_j_it){
               string x_j = *x_j_it;
               if (x_i != x_j){
                  if (p[l_p][x_j]!=x_i){
                     g[x_i][x_j] = f[l_p][x_j] + distance_dict[x_j][x_i];
                     g_type[x_i][x_j] = q_route[l_p][x_j];
                     g_type[x_i][x_j].push_back(x_i);
                  } else {
                     g[x_i][x_j] = phi[l_p][x_j] + distance_dict[x_j][x_i];
                     g_type[x_i][x_j] = q_route_2[l_p][x_j];
                     g_type[x_i][x_j].push_back(x_i);
                  }
               }
            }
         }
      }
      for (vector<string>::iterator x_i_it = N.begin(); x_i_it != N.end(); ++x_i_it){
         string x_i = *x_i_it;
         int q_p = Q - quantities[x_i];
         if (q_p > 0){
            // Find the minimum
            string arg_min_1;
            string arg_min_2;
            double min_1 = inf;
            double min_2 = inf;
            for (vector<string>::iterator x_j_it = N.begin(); x_j_it != N.end(); ++x_j_it){
               string x_j = *x_j_it;
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
            q_route[l][x_i] = g_type[x_i][arg_min_1];
            q_route_2[l][x_i] = g_type[x_i][arg_min_2];
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
