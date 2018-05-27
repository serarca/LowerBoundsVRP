#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <limits>
#include <iterator>
#include <numeric>
#include "lower_bounds.h"
#include <algorithm>
#include <numeric>


using namespace std;

PossibleValues possible_values(vector<int>& quantities, int truck_capacity){
   vector<int> values(truck_capacity);

   for (int i=0; i< truck_capacity; i++){
      values[i] = i + 1;
   }
   map<int, int> values_pos;
   for (int i = 0; i < (int) values.size(); i++){
      values_pos[values[i]] = i;
   }
   PossibleValues possible;
   possible.values = values;
   possible.values_pos = values_pos;
   return possible;
}

// Write lower bounds
LowerBound lower_bound_(
   vector<int> H,
   vector<int> capacities,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> distance_dict,
   vector<double> mu,
   vector<double> lamb
){

   //Calculate lengths
   int len_N = N.size();
   int len_H = H.size();


   //The vectors of minima
   vector<vector<double>> b(len_N, vector<double>(len_H));
   vector<vector<int>> val(len_N, vector<int>(len_H));
   vector<vector<vector<int>>> b_routes(len_N, vector<vector<int>>(len_H, vector<int>(0)));

   for (int h = 0; h < len_H; h++){
      int truck_capacity = capacities[h];
      PossibleValues possible = possible_values(quantities, truck_capacity);

      QRoutes qroutes = construct_q_routes_(H[h], truck_capacity, N, distance_dict, possible.values, possible.values_pos, quantities);


      // We find the minimum l

      int len_values = possible.values.size();
      for (int n = 0; n < len_N; n++){
         vector<double> b_values(len_values);

         for (int l = 0; l < len_values; l++){
            b_values[l] = (qroutes.psi[l][n] - mu[h]) * (double)quantities[n]/(double)possible.values[l];
         }

         std::vector<double>::iterator l_it = std::min_element(b_values.begin(), b_values.end());
         int l_min = std::distance(b_values.begin(), l_it);
         b[n][h] = b_values[l_min];
         b_routes[n][h] = qroutes.psi_route[l_min][n];
         val[n][h] = possible.values[l_min];


      }


   }

   // The second vector of minima
   vector<double> b_min(len_N);
   vector<int> val_min(len_N);
   vector<int> h_min(len_N);
   vector<vector<int>> b_min_routes(len_N, vector<int>(0));
   for (int n = 0; n < len_N; n++){
      std::vector<double>::iterator h_it = std::min_element(b[n].begin(), b[n].end());
      int h_m = std::distance(b[n].begin(), h_it);
      b_min[n] = b[n][h_m];
      b_min_routes[n] = b_routes[n][h_m];
      val_min[n] = val[n][h_m];
      h_min[n] = h_m;


   }

   //Calculate number of visits to each node
   vector<vector<int>> visits(len_N,vector<int>(len_N,0));
   for (int n = 0; n < len_N; n++){
      for(vector<int>::iterator it_route = b_min_routes[n].begin();  it_route!=b_min_routes[n].end(); ++it_route){
         if (*it_route < len_N){
            visits[n][*it_route] += 1;
         }
      }
   }


   // Construct theta
   vector<double> theta(len_N,0);
   for (int j = 0; j<len_N; j++){
      for (int i=0; i<len_N; i++){
         theta[j] += ((double) quantities[i])/((double) val_min[i])*(double)visits[i][j];
      }
   }
   // Construct rho
   vector<double> rho(len_H,0);
   for (int n=0; n<len_N; n++){
      rho[h_min[n]] += ((double) quantities[n])/((double) val_min[n]);
   }

   // Construct dual variables
   vector<double> u(len_N,0);
   for (int n=0; n<len_N; n++){
      u[n] = b_min[n] + lamb[n];
   }
   double start = 0;

   double z_lb = std::accumulate(u.begin(), u.end(), start) + std::accumulate(mu.begin(), mu.end(), start);

   LowerBound lb;
   lb.z_lb = z_lb;
   lb.u = u;
   lb.theta = theta;
   lb.rho = rho;

   return lb;

}


// This function returns the q-paths that the q-routes function uses
// to calulate lower bounds
QPaths construct_q_paths_(
   int h,
   int truck_capacity,
   vector<int> N,
   vector<vector<double>> distance_dict,
   vector<int> values,
   map<int,int> values_pos,
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


QRoutes construct_q_routes_(
   int h,
   int truck_capacity,
   vector<int> N,
   vector<vector<double>> distance_dict,
   vector<int> values,
   map<int,int> values_pos,
   vector<int> quantities
){

   QPaths qpaths_l = construct_q_paths_(h,truck_capacity,N,distance_dict,values,values_pos,quantities,"left");
   QPaths qpaths_r = construct_q_paths_(h,truck_capacity,N,distance_dict,values,values_pos,quantities,"right");

   vector<vector<double>>& f_l = qpaths_l.f;
   vector<vector<double>>& phi_l = qpaths_l.phi;
   vector<vector<int>>& p_l = qpaths_l.p;
   vector<vector<vector<int>>>& q_route_l = qpaths_l.q_route;
   vector<vector<vector<int>>>& q_route_2_l = qpaths_l.q_route_2;

   vector<vector<double>>& f_r = qpaths_r.f;
   vector<vector<double>>& phi_r = qpaths_r.phi;
   vector<vector<int>>& p_r = qpaths_r.p;
   vector<vector<vector<int>>>& q_route_r = qpaths_r.q_route;
   vector<vector<vector<int>>>& q_route_2_r= qpaths_r.q_route_2;




   //Calculate lengths
   int len_values = values.size();
   int len_N = N.size();

   //Construct infinity
   double inf = numeric_limits<double>::infinity();

   vector<vector<double>> psi(len_values, vector<double> (len_N));
   vector<vector<vector<int>>> psi_route(len_values, vector<vector<int>> (len_N, vector<int>(0)));

   for(int l = 0; l < len_values; l++) {
      for(int n = 0; n < len_N; n++) {
         int min_w = quantities[n];
         int max_w = values[l] + quantities[n];
         double min_val = inf;
         vector<int> min_coord(3,-1);
         double val(0);
         int option;
         for (int l_1 = 0; l_1 < len_values; l_1++){
            int q = values[l_1];
            if (q>=min_w && q<max_w){
               int l_2 = values_pos[values[l] + quantities[n] - q];
               if (p_l[l_1][n]!=p_r[l_2][n] || (p_l[l_1][n] == p_r[l_2][n] && p_l[l_1][n]==h)){
                  val = f_l[l_1][n]+f_r[l_2][n];
                  option = 0;
               } else{
                  if (f_l[l_1][n]+phi_r[l_2][n]<phi_l[l_1][n]+f_r[l_2][n]){
                     val = f_l[l_1][n]+phi_r[l_2][n];
                     option = 1;
                  } else {
                     val = phi_l[l_1][n]+f_r[l_2][n];
                     option = 2;
                  }
               }
               if (val < min_val){
                  min_val = val;
                  min_coord[0] = option;
                  min_coord[1] = l_1;
                  min_coord[2] = l_2;
               }
            }
         }
         psi[l][n] = min_val;
         if (min_coord[0]!=-1){
            if (min_coord[0] == 0){
               vector<int> & A = q_route_l[min_coord[1]][n];
               vector<int> & B = q_route_r[min_coord[2]][n];
               psi_route[l][n].reserve( A.size() + B.size() - 1 ); // preallocate memory
               psi_route[l][n].insert( psi_route[l][n].end(), A.begin(), A.end() );
               psi_route[l][n].insert( psi_route[l][n].end(), B.rbegin() + 1, B.rend() ); // insert backwards
            } else if (min_coord[0] == 1){
               vector<int> & A = q_route_l[min_coord[1]][n];
               vector<int> & B = q_route_2_r[min_coord[2]][n];
               psi_route[l][n].reserve( A.size() + B.size() - 1 ); // preallocate memory
               psi_route[l][n].insert( psi_route[l][n].end(), A.begin(), A.end() );
               psi_route[l][n].insert( psi_route[l][n].end(), B.rbegin() + 1, B.rend() ); // insert backwards
            } else if (min_coord[0] == 2){
               vector<int> & A = q_route_2_l[min_coord[1]][n];
               vector<int> & B = q_route_r[min_coord[2]][n];
               psi_route[l][n].reserve( A.size() + B.size() - 1 ); // preallocate memory
               psi_route[l][n].insert( psi_route[l][n].end(), A.begin(), A.end() );
               psi_route[l][n].insert( psi_route[l][n].end(), B.rbegin() + 1, B.rend() ); // insert backwards
            }
         }
      }
   }

   QRoutes qroutes;
   qroutes.psi = psi;
   qroutes.psi_route = psi_route;

   return qroutes;
}
