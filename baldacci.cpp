#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <limits>
#include <iterator>
#include <numeric>
#include "baldacci.h"
#include <algorithm>
#include <numeric>
#include <list>
#include <tuple>

using GenPath = vector<list<Path>>;
using GenRoute = vector<list<Route>>;

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
                  T[i].erase(p);
               } else {
                  ++p;
               }
            }
         }
      }
   }

   return P;
}

GenRoute GENROUTE(
   int Delta,
   double gamma,
   int h,
   int capacity,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> distance_dict){

   int len_N = N.size();

   GenPath P_l = GENPATH(Delta, gamma, h, capacity, N, quantities, distance_dict, "left");
   GenPath P_r = GENPATH(Delta, gamma, h, capacity, N, quantities, distance_dict, "right");

   vector<list<Route>> T(len_N, list<Route>(0));
   vector<list<Route>> R(len_N, list<Route>(0));

   // Set of added pairs of routes
   vector<set<tuple<int, int>>> added(len_N);

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
      if (iterations>0)
         break;
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


      bool valid = true;
      if (total_load > capacity){
         // Check if violate capacity
         valid = false;
      } else if (((double) load_l < (double)total_load/2.0) or
                  ((double) load_r < (double)total_load/2.0) or
                  ((double) load_l > (double)total_load/2.0 +(double) quantities[arg_min]) or
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
         }
      }
      // Order by cost
      std::sort (new_routes.begin(), new_routes.end(), compare_routes);
      // Insert them
      int loc = 0;
      auto p = T[arg_min].begin();
      while (p!=T[arg_min].end()){
         if (loc == 2)
            break;
         if ((p->cost)>new_routes[loc].cost){
            T[arg_min].insert(p,new_routes[loc]);
            loc++;
         } else {
            ++p;
         }
      }
      // Finish inserting
      while (loc != 2){
         T[arg_min].insert(p,new_routes[loc]);
         loc++;
      }

   }
   return R;

}
