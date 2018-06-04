#pragma once


#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <limits>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <list>
#include "lower_bounds.h"
#include "baldacci.h"

#include<cmath>
#include "prettyprint.hpp"




using GenPath = vector<list<Path>>;
using GenRoute = vector<list<SimpleRoute>>;

bool compare_routes (Route i, Route j) { return (i.cost<j.cost); }

vector<list<Path>> GENPATH(
   int Delta,
   double gamma,
   int h,
   int capacity,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> distance_dict,
   string direction){

   int len_N = N.size();

   vector<list<Path>> P(N.size() + 1, list<Path>(0));
   vector<list<Path>> T(N.size() + 1, list<Path>(0));

   Path init;
   init.path.push_front(h);
   init.cost = 0;
   init.lower_bound = 0;
   init.load = 0;
   init.end = h;
   T[N.size()].push_front(init);

   // Define infinity
   double inf = numeric_limits<double>::infinity();
   int count_paths = 0;
   while (true){
      // Check the route with the smallest cost
      int arg_min = -1;
      double cost_min = inf;
      // No typo here with the <=
      for(int i = 0; i <= len_N; i++){
         if (T[i].size() > 0){
            double cost = T[i].front().cost;
            if (cost<cost_min){
               cost_min = cost;
               arg_min = i;
            }
         }
      }
      // Break if no more routes
      if (arg_min == -1)
         break;
      // Check the first route and add it to P
      P[arg_min].splice(P[arg_min].end(),T[arg_min],T[arg_min].begin());
      count_paths += 1;
      // Break if too many paths
      if (count_paths==Delta)
         break;
      // Extract the element
      Path p_star = P[arg_min].back();
      // If path violates capacity, go to the next one
      if ((double) p_star.load > ((double) capacity)/2.0)
         continue;
      // Add new paths
      for (int i = 0; i < len_N; i++){
         // Check node is not in path
         auto it = p_star.nodes.find(N[i]);
         if (it == p_star.nodes.end()){
            Path new_path;
            new_path.path = p_star.path;
            new_path.path.push_back(i);
            new_path.nodes = p_star.nodes;
            new_path.nodes.insert(i);
            new_path.load = p_star.load + quantities[i];
            new_path.end = N[i];
            if (direction == "left"){
               new_path.cost = p_star.cost + distance_dict[p_star.end][i];
               new_path.lower_bound = new_path.cost;
            }
            if (direction == "right"){
               new_path.cost = p_star.cost + distance_dict[i][p_star.end];
               new_path.lower_bound = new_path.cost;
            }
            //Check if the new path has a cost too high
            if (new_path.lower_bound >= gamma)
               continue;
            //Check if the new path has a load too high
            if (new_path.load > capacity)
               continue;
            //Check if this new path is dominated by any path in P
            bool dominated = false;
            for (auto p = P[i].begin(); p!=P[i].end(); ++p){
                if ((p->end == new_path.end) && (p->cost <= new_path.cost) && (p->nodes == new_path.nodes)){
                   dominated = true;
                   break;
                }
            }
            if (dominated)
                continue;

            // Check if the path is dominated by any path in T
            auto p = T[i].begin();
            for ( ; p!=T[i].end(); ++p){
               if ((p->end == new_path.end) && (p->cost <= new_path.cost) && (p->nodes == new_path.nodes)){
                  dominated = true;
                  break;
               }
               if (p->cost > new_path.cost)
                  break;
            }
            if (dominated)
               continue;
            // Insert the path
            T[i].insert(p,new_path);
            // Delete dominated elements
            while (p!=T[i].end()){
               if ((p->end == new_path.end) && (p->cost > new_path.cost) && (p->nodes == new_path.nodes)){
                  p=T[i].erase(p);
               } else {
                  ++p;
               }
            }
         }
      }
   }

   return P;
}

list<SimpleRoute> GENROUTE(
   int Delta,
   double gamma,
   int h,
   int capacity,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> distance_dict,
   vector<vector<double>> geo_distance){

   int len_N = N.size();

   GenPath P_l = GENPATH(Delta, gamma, h, capacity, N, quantities, distance_dict, "left");
   GenPath P_r = GENPATH(Delta, gamma, h, capacity, N, quantities, distance_dict, "right");

   vector<list<Route>> T(len_N, list<Route>(0));
   vector<list<Route>> R(len_N, list<Route>(0));

   // Set of added pairs of routes
   vector<set<std::tuple<int, int>>> added(len_N);

   // Initialize the routes
   for (int i=0; i<len_N; i++){
      added[i].insert(std::make_tuple(-1, -1));
      if (P_l[i].size()>0 && P_r[i].size()>0){
         Route init;
         init.path_l = P_l[i].begin();
         init.path_r = P_r[i].begin();
         init.cost = (init.path_l->cost)+ (init.path_r->cost);
         init.index_l = 0;
         init.index_r = 0;
         T[i].push_front(init);
         added[i].insert(std::make_tuple(0, 0));
      }
   }
   // Define infinity
   double inf = numeric_limits<double>::infinity();
   int iterations = 0;
   while(true){
      ++iterations;
      // Check the route with the smallest cost
      int arg_min = -1;
      double cost_min = inf;
      for(int i = 0; i < len_N; i++){
         if (T[i].size() > 0){
            double cost = T[i].front().cost;
            if (cost<cost_min){
               cost_min = cost;
               arg_min = i;
            }
         }
      }
      // Break if no more routes
      if (arg_min == -1)
         break;
      // Break if exceed max cost
      if (cost_min>gamma){
         break;
      }
      // Check the first route and pop it
      Route r_star = T[arg_min].front();
      T[arg_min].pop_front();

      // Fill the load
      int load_l = r_star.path_l->load;
      int load_r = r_star.path_r->load;
      int total_load = load_l + load_r - quantities[arg_min];
      r_star.load = total_load;

      //Fill the median
      r_star.median = arg_min;
      bool valid = true;
      if (total_load > capacity){
         // Check if violate capacity
         valid = false;
      } else if (((double) load_l < (double)total_load/2.0) ||
                  ((double) load_r < (double)total_load/2.0) ||
                  ((double) load_l > (double)total_load/2.0 +(double) quantities[arg_min]) ||
                  ((double) load_r > (double)total_load/2.0 +(double) quantities[arg_min])){
         // Check if violate median
         valid = false;
      } else {
         // Check if violate empty intersection
         set<int> intersection;
         set<int> s1 = r_star.path_l->nodes;
         set<int> s2 = r_star.path_r->nodes;
         set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(), std::inserter(intersection,intersection.begin()));
         int myints[]= {arg_min};
         set<int> comparison (myints,myints+1);
         if (intersection != comparison){
            valid = false;
         }
      }
      if (valid){
         // Calculate the union of the nodes in the path
         set<int> s1 = r_star.path_l->nodes;
         set<int> s2 = r_star.path_r->nodes;
         set_union(s1.begin(),s1.end(),s2.begin(),s2.end(), std::inserter(r_star.nodes,r_star.nodes.begin()));

         for (int i = 0; i<len_N; i++){
            for (auto r = R[i].begin(); r!=R[i].end(); ++r){
               if ((r->nodes) == (r_star.nodes)){
                  valid = false;
                  break;
               }
            }
            if (!valid)
               break;
         }
      }
      // Push if valid
      if (valid){
         R[arg_min].push_back(r_star);
      }
      // We add new routes
      vector<Route> new_routes;
      if (r_star.index_l + 1 < (int)P_l[arg_min].size()){
         // Check if not added already
         if (added[arg_min].count(std::make_tuple(r_star.index_l + 1, r_star.index_r)) == 0) {
            Route n_route;
            n_route.path_l = std::next(r_star.path_l,1);
            n_route.path_r = r_star.path_r;
            n_route.cost = (n_route.path_l->cost)+ (n_route.path_r->cost);
            n_route.index_l = r_star.index_l + 1;
            n_route.index_r = r_star.index_r;
            new_routes.push_back(n_route);
	    added[arg_min].insert(std::make_tuple(n_route.index_l, n_route.index_r));
         }
      }
      if (r_star.index_r + 1 < (int)P_r[arg_min].size()){
         if (added[arg_min].count(std::make_tuple(r_star.index_l, r_star.index_r + 1)) == 0){
            Route n_route;
            n_route.path_r = std::next(r_star.path_r,1);
            n_route.path_l = r_star.path_l;
            n_route.cost = (n_route.path_l->cost)+ (n_route.path_r->cost);
            n_route.index_r = r_star.index_r + 1;
            n_route.index_l = r_star.index_l;
            new_routes.push_back(n_route);
            added[arg_min].insert(std::make_tuple(n_route.index_l, n_route.index_r));
         }
      }
      // Order by cost
      std::sort (new_routes.begin(), new_routes.end(), compare_routes);
      // Insert them
      int loc = 0;
      auto p = T[arg_min].begin();
      while (p!=T[arg_min].end()){
         if (loc == (int) new_routes.size())
            break;
         if ((p->cost)>new_routes[loc].cost){
            T[arg_min].insert(p,new_routes[loc]);
            loc++;
         } else {
            ++p;
         }
      }
      // Finish inserting
      while (loc != (int) new_routes.size()){
         T[arg_min].insert(p,new_routes[loc]);
         loc++;
      }
   }


   // Take the routes and reconstruct them, since we do not want pointers to
   // iterators but the paths themselves
   std::list<SimpleRoute> SimpleRoutes;

   for (int i = 0; i < len_N; i++){
      for (auto route:R[i]){
         SimpleRoute copy;
         // Add first path
         for (auto it = (route.path_l->path).begin(); it!=(route.path_l->path).end(); it++){
            copy.path.push_back(*it);
         }
         // Add second path in opposite order (make sure not to add first element)
         for (auto it = (route.path_r->path).rbegin(); it!=(route.path_r->path).rend(); it++){
            if (it!=(route.path_r->path).rbegin())
               copy.path.push_back(*it);
         }
         copy.index_l = route.index_l;
         copy.index_r = route.index_r;
         copy.cost = route.cost;
         copy.load = route.load;
         copy.median = route.median;
         copy.truck = h;
         // Fill the geo_cost
         copy.geo_cost = 0;
         auto it1 = copy.path.begin();
         auto it2 = copy.path.begin();
         it2++;
         while(it2!=copy.path.end()){
            copy.geo_cost+=geo_distance[*it1][*it2];
            it1++;
            it2++;
         }
         SimpleRoutes.push_back(copy);
      }
   }

   for (int i = 0; i < len_N; i++){
      if (R[i].size() == 0){
         SimpleRoute n_route;
         n_route.path.push_back(h);
         n_route.path.push_back(i);
         n_route.path.push_back(h);
         n_route.cost = distance_dict[i][h]+ distance_dict[h][i];
         n_route.index_r = -1;
         n_route.index_l = -1;
         n_route.load = quantities[i];
         n_route.median = i;
         n_route.truck = h;
         // Fill the geo_cost
         n_route.geo_cost = 0;
         auto it1 = n_route.path.begin();
         auto it2 = n_route.path.begin();
         it2++;
         while(it2!=n_route.path.end()){
            n_route.geo_cost+=geo_distance[*it1][*it2];
            it1++;
            it2++;
         }
         SimpleRoutes.push_back(n_route);
      }
   }

   // Calculate the geo_costs


   return SimpleRoutes;

}

// This function takes a set of routes with given reduced costs and finds
// the lower bound they induce, as well as updates the reduced costs and
// returns a new version of mu and lambda
LowerBound lower_bound_M2(
   vector<list<SimpleRoute>> Routes,
   vector<int> H,
   vector<int> capacities,
   vector<int> N,
   vector<int> quantities,
   vector<double> mu,
   vector<double> lamb
){


   //Calculate lengths
   int len_N = N.size();
   int len_H = H.size();

   // Infinity
   double inf = numeric_limits<double>::infinity();

   // Update the costs of the routes
   for (auto & truck_routes:Routes){
      for (auto & route:truck_routes){
         route.cost = route.geo_cost;
         route.cost -= mu[route.truck - len_N];
         for (auto node:route.path){
            if (node<len_N){
               route.cost -= lamb[node];
            }
         }
      }
   }
   /*
   for (auto & truck_routes:Routes){
      for (auto & route:truck_routes){
         cout<<"Cost "<<route.cost<<" Geo "<<route.geo_cost<<endl;
      }
   }
*/
   //The vectors of minima
   vector<vector<double>> b(len_N, vector<double>(len_H, inf));
   vector<vector<SimpleRoute>> b_routes(len_N, vector<SimpleRoute>(len_H));

   for (int h = 0; h < len_H; h++){
      for (auto route:Routes[h]){
         for (auto node:route.path){
            // Verify the node is a farmer
            if (node<len_N){
               double new_b = route.cost / ((double) route.load) * ((double)quantities[node]);
               if (new_b < b[node][h]){
                  b[node][h] = new_b;
                  b_routes[node][h] = route;
               }
            }
         }
      }

   }


   // The second vector of minima
   vector<double> b_min(len_N);
   vector<SimpleRoute> b_min_routes(len_N);
   for (int n = 0; n < len_N; n++){
      std::vector<double>::iterator h_it = std::min_element(b[n].begin(), b[n].end());
      int h_m = std::distance(b[n].begin(), h_it);
      b_min[n] = b[n][h_m];
      b_min_routes[n] = b_routes[n][h_m];
   }

   //Calculate number of visits to each node
   vector<vector<int>> visits(len_N,vector<int>(len_N,0));
   for (int n = 0; n < len_N; n++){
      for(auto it_route = b_min_routes[n].path.begin();  it_route!=b_min_routes[n].path.end(); ++it_route){
         if (*it_route < len_N){
            visits[n][*it_route] += 1;
         }
      }
   }


   // Construct theta
   vector<double> theta(len_N,0);
   for (int j = 0; j<len_N; j++){
      for (int i=0; i<len_N; i++){
         theta[j] += ((double) quantities[i])/((double) b_min_routes[i].load)*(double)visits[i][j];
      }
   }


   // Construct rho
   vector<double> rho(len_H,0);
   for (int n=0; n<len_N; n++){
      rho[b_min_routes[n].truck - len_N] += ((double) quantities[n])/((double) b_min_routes[n].load);
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

// Given a set of routes, calculate the best lower bounds that they generate
// We initialize at mu and lamb
DualSolution lower_bound_optimizer_M2(
   int iterations,
   vector<list<SimpleRoute>> Routes,
   double z_ub,
   double epsilon,
   vector<int> H,
   vector<int> capacities,
   vector<int> N,
   vector<int> quantities,
   vector<double> mu,
   vector<double> lamb){

   //Calculate lengths
   int len_N = N.size();
   int len_H = H.size();
   // Define infinity
   double infinity = numeric_limits<double>::infinity();
   // Define u
   vector<double> u(len_N);

   // Vectors to store the optimal values
   vector<double> lamb_opt(len_N);
   vector<double> u_opt(len_N);
   vector<double> v_opt(len_H);
   // Here we store the values of the iterations
   vector<double> values;
   double max_val = -infinity;
   vector<double> save_theta;
   vector<double> save_rho;
   double g_den_1 = 0;
   double g_den_2 = 0;
   double g = 0;
   for (int iteration = 0; iteration<iterations; iteration++){
      // We pass these routes to the algorithm that calculates the lower bound
      LowerBound lb = lower_bound_M2(Routes, H, capacities, N, quantities, mu, lamb);

      // Check if the lower bound that we get improves
      if (lb.z_lb > max_val){
         max_val = lb.z_lb;
         u_opt = lb.u;
         v_opt = mu;
         lamb_opt = lamb;
      }
      cout<<lb.z_lb<<endl;
      //cout<<lamb<<endl;
      //cout<<mu<<endl;
      save_theta = lb.theta;
      save_rho = lb.rho;
      values.push_back(lb.z_lb);
      // We calculate g for the step of the algorithm
      g_den_1 = 0;
      g_den_2 = 0;
      for (int i = 0; i< len_N; i++){
         g_den_1 += pow(lb.theta[i]-1.0,2);
      }
      for (int i = 0; i< len_H; i++){
         g_den_2 += pow(lb.rho[i]-1.0,2);
      }
      g = (z_ub - lb.z_lb)/(g_den_1 + g_den_2);
      // Update lambda and mu
      for (int i = 0; i<len_N; i++){
         lamb[i] = lamb[i] - epsilon*g*(lb.theta[i]-1.0);
      }
      for (int i = 0; i<len_H; i++){
         mu[i] = min(mu[i] - epsilon*g*(lb.rho[i]-1.0),0.0);
      }

      // Check that we are not getting exploting reduced variables
      double explotion = 0;
      for (int i = 0; i< len_N; i++){
         explotion+=fabs(lamb[i]);
      }
      if (explotion > pow(10,16)){
         cout<<"Diverging reduced variables"<<endl;
         break;
      }

      // We update epsilon
      if ((int) values.size() >= 7){
         // Calculate changes
         vector<double> grad;
         for (int i = 0; i < (int)values.size()-1; i++){
            grad.push_back(values[i+1] - values[i]);
         }
         // Calculate jumps
         vector<int> jumps;
         for (int i = 0; i < (int)grad.size()-1; i++){
            jumps.push_back(signbit(grad[i+1])!=signbit(grad[i]));
         }
         // If too many jumps reduce epsilon
         int n_jumps = 0;
         for (int i = (int)jumps.size()-5; i<(int)jumps.size(); i++ ){
            n_jumps += jumps[i];
         }
         if (n_jumps >= 3){
            epsilon = epsilon/1.5;
            cout<<"New epsilon "<<epsilon<<endl;
            values.clear();
         }
      }
      // Check if we have reached a zero gradient
      double gradient_norm = 0;
      for (int i = 0; i < len_N; i++){
         gradient_norm += pow(lb.theta[i] - 1.0,2);
      }
      for (int i = 0; i < len_H; i++){
         gradient_norm += pow(lb.rho[i] - 1.0,2);
      }

      if (gradient_norm < pow(10.0, -20)){
         cout<<"Reached zero gradient"<<endl;
         break;
      }
   }

   DualSolution new_bound;
   new_bound.z_lb = max_val;
   new_bound.u = u_opt;
   new_bound.v = v_opt;
   new_bound.lamb = lamb_opt;
   cout<<"New opt "<< max_val<<endl;
   return new_bound;
}


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
){
   //Calculate lengths
   int len_N = N.size();
   int len_H = H.size();
   // Define infinity
   double infinity = numeric_limits<double>::infinity();
   // Define u
   vector<double> u(len_N);

   // Define the vector of routes for each truck
   vector<list<SimpleRoute>> Routes(len_H);

   // Calculate reduced costs
   vector<vector<double>> distance_dict = reduced_cost_matrix(geo_distance, lamb, mu);

   // We start by generating routes for all of the trucks
   for (auto h:H){
      Routes[h - len_N] = GENROUTE(Delta, gamma, h, capacities[h - len_N], N, quantities, distance_dict, geo_distance);
   }
   // Vectors to store the optimal values
   vector<double> lamb_opt(len_N);
   vector<double> u_opt(len_N);
   vector<double> v_opt(len_H);
   // Here we store the values of the iterations
   vector<double> values;

   while(true){
      // We pass these routes to the algorithm that calculates the lower bound
      DualSolution ds = lower_bound_optimizer_M2(sub_iterations, Routes, z_ub, epsilon, H, capacities, N, quantities, mu, lamb);
      vector<list<SimpleRoute>> newRoutes(len_H);
      vector<vector<double>> new_distance_dict = reduced_cost_matrix(geo_distance, ds.lamb, ds.v);
      bool empty = true;
      for (auto h:H){
         newRoutes[h - len_N] = GENROUTE(Delta_zero, gamma_zero, h, capacities[h - len_N], N, quantities, new_distance_dict, geo_distance);
         for (auto route:newRoutes[h - len_N]){
            if (route.cost < 0){
               Routes[h - len_N].push_back(route);
               empty = false;
            }
         }
      }
      if (empty){
         break;
      }

   }
}
