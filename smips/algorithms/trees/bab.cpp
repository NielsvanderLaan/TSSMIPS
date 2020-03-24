#include "tree.h"

vector<double> Tree::bab(bool affine, double tol)
{ 
  vector<double> incumbent(d_problem.d_n1);
  double local_tol = tol;      // conditional on integer first-stage
  int max_iter = 0;            // limit on number of consecutive Benders iterations without integer first-stage solutions
  
  while (d_UB_global > d_LB_global + tol)
  {   
    size_t node_idx = distance(d_LB_nodes.begin(), min_element(d_LB_nodes.begin(), d_LB_nodes.end()));
    cout << "\nExploring node " << node_idx << '\n';

    bool fathom = solve(node_idx, incumbent, affine, tol, local_tol, max_iter);  // solve() also updates global bounds
        
    cout << "GLOBAL LB = " << d_LB_global << " GLOBAL UB = " << d_UB_global << '\n';
    
    local_tol = max(tol / 10, local_tol / 1.1); 

      
    if (fathom)          // node is infeasible or LB(node) > global_UB
      continue;          // else: branch
      
    Split split = branch_var(node_idx);
    if (split.var == -1)
    {
      if (d_LB_nodes[node_idx] <= d_LB_global)  // no improvement possible
        break;
      else                                      // global improvement possible  
        continue;
    }  
    
    branch(node_idx, split);
  }
  
  if (d_UB_global > d_LB_global + tol)
    cout << "GAP > tolerance. GAP = " << d_UB_global - d_LB_global << '\n';
  
  
  cout << "number of nodes: " << d_nodes.size() << '\n';
  return incumbent;
}

    
    






