#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <limits>
#include <iterator>

using namespace std;
// The struct of results for q-routes
struct QRoutes {
   vector<vector<double>> psi;
   vector<vector<vector<int>>> psi_route;
};

// The struct of results
struct QPaths {
   vector<vector<double>> f;
   vector<vector<double>> phi;
   vector<vector<int>> p;
   vector<vector<vector<int>>> q_route;
   vector<vector<vector<int>>> q_route_2;
};




// This function returns the q-paths that the q-routes function uses
// to calulate lower bounds
QPaths construct_q_paths_(
   int& h,
   int& truck_capacity,
   vector<int>& N,
   vector<vector<double>>& distance_dict,
   vector<int>& values,
   map<double,int>& values_pos,
   vector<int>& quantities,
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
   map<double,int> values_pos,
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
