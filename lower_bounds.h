#include <vector>
#include <map>
#include <string>
#include <iostream>


using namespace std;


// The struct of results
struct QPaths {
   map<int,map<string,double>> f;
   map<int,map<string,double>> phi;
   map<int,map<string,string>> p;
   map<int,map<string,vector<string>>> q_route;
   map<int,map<string,vector<string>>> q_route_2;
};

QPaths construct_q_paths_(string h,
   int truck_capacity,
   vector<string> N,
   map<string,map<string,double>> distance,
   vector<int> values,
   map<double,int> values_pos,
   map<string,int> quantities,
   string direction);
