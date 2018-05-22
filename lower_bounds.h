#include <vector>
#include <map>
#include <string>
#include <iostream>

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

// The struct of lower bounds
struct LowerBound {
   double z_lb;
   vector<double> theta;
   vector<double> rho;
   vector<double> u;
};

// A struct with posible values and their inverse maps
struct PossibleValues {
   vector<int> values;
   map<int,int> values_pos;
};

PossibleValues possible_values(vector<int>& quantities, int truck_capacity);

LowerBound lower_bound_(
   vector<int> H,
   vector<int> capacities,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> distance_dict,
   vector<double> mu,
   vector<double> lamb
);

QPaths construct_q_paths_(
   int h,
   int truck_capacity,
   vector<int> N,
   vector<vector<double>> distance_dict,
   vector<int> values,
   map<int,int> values_pos,
   vector<int> quantities,
   string direction
);

QRoutes construct_q_routes_(
   int h,
   int truck_capacity,
   vector<int> N,
   vector<vector<double>> distance_dict,
   vector<int> values,
   map<int,int> values_pos,
   vector<int> quantities
);
