#include <vector>
#include <map>
#include <string>
#include <iostream>


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
);
