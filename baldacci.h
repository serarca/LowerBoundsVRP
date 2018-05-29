#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <list>
#include <set>

using namespace std;

//The struct of paths
struct Path{
   list<int> path;
   set<int> nodes;
   double cost;
   double lower_bound;
   int load;
   int end;
};

//The struct of paths
struct Route{
   list<Path>::iterator path_l;
   list<Path>::iterator path_r;
   int index_l;
   int index_r;
   set<int> nodes;
   double cost;
   int load;
   int median;
};
