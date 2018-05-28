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
      Path p_star = P[arg_min].front();
      // If path violates capacity, go to the next one
      cout<<p_star.load<<endl;
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

   return T;
}
