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
