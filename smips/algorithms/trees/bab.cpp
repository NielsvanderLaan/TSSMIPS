#include "tree.h"

vector<double> Tree::bab(bool lp_cuts, bool sb_cuts, bool zk_cuts, bool strong_cuts, bool affine, double tol)
{ 
  vector<double> incumbent(d_problem.d_n1);
  double local_tol = tol / 10;

  size_t tree_size = 1;
  while (d_UB_global > d_LB_global + tol && not d_nodes.empty())
  {
    cout << "\nLBs: ";
    for_each(d_LB_nodes.begin(), d_LB_nodes.end(), [](double val){ cout << val << ' '; });
    size_t node_idx = distance(d_LB_nodes.begin(), min_element(d_LB_nodes.begin(), d_LB_nodes.end()));
    cout << "\nExploring node " << node_idx << endl;

    bool branch = solve(lp_cuts, sb_cuts, zk_cuts, strong_cuts, affine, node_idx, incumbent, local_tol);  // solve() also updates global bounds
    cout << "GLOBAL LB = " << d_LB_global << " GLOBAL UB = " << d_UB_global << endl;
    
    //local_tol = max(tol / 10, local_tol / 1.1);
    fathom();

    if (not branch) continue;   // else branch
      
    Split split = branch_var(node_idx);
    if (split.var == -1)
    {
      if (d_LB_nodes[node_idx] <= d_LB_global + tol)  // no improvement possible
        break;
      else                                      // global improvement possible  
        continue;
    }
    add_branch(node_idx, split);
    ++tree_size;
  }
  
  //if (d_UB_global > d_LB_global + tol)
  //  cout << "GAP > tolerance. GAP = " << d_UB_global - d_LB_global << '\n';

  cout << "global LB: " << d_LB_global << '\n';
  cout << "global UB: " << d_UB_global << '\n';
  cout << "number of nodes: " << d_nodes.size() << '\n';
  cout << "tree size: " << tree_size << '\n';
  return incumbent;
}

    
    






