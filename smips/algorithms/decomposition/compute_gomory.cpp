#include "benders.h"

double Benders::compute_gomory(size_t s, vector<int> &vbasis, vector<int> &cbasis, vector<double> &ws,vector<double> &alpha)
{
  vector<double> basis(d_n2 + d_m2);
  copy(vbasis.begin(), vbasis.begin() + d_n2, basis.begin());
  copy(cbasis.begin(), cbasis.begin() + d_m2, basis.begin() + d_n2);
  
  vector<vector<double>> &visited_bases = d_visited[s];
  
  auto end = visited_bases.end();
  auto begin = visited_bases.begin();
  auto it = find(begin, end, basis);

  double gom_obj;
  if (it == end)  // not visited, update and solve d_gomory
  {
    vector<double> omega_alpha(d_m2);                       // rhs vector of gomory relaxation
    transform(omega_alpha.begin(), omega_alpha.end(), alpha.begin(), omega_alpha.end(), minus<double>());

    d_gomory.update(omega_alpha, vbasis, cbasis);  // update gomory relaxation
  
    gom_obj = d_gomory.solve();                    // solve gomory relaxation and store objective bound 

    visited_bases.push_back(basis);                // basis is now visited: store it 
    d_objectives[s].push_back(gom_obj);            // store corresponding objective value   
  } else          // visited before
    gom_obj = d_objectives[s][distance(begin, it)];                // and retrieve corresponding objective value

  return gom_obj;
}
