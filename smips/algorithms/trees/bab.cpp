#include "tree.h"

vector<double> Tree::bab(vector<Type> types, bool rcuts, bool fenchel, size_t max_rounds, double tol, double time_limit)
{
  vector<double> incumbent(d_problem.d_n1);

  size_t tree_size = 1;
  auto t1 = chrono::high_resolution_clock::now();
  while (d_UB_global > d_LB_global + tol && not d_nodes.empty())
  {
    double elapsed = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t1).count() / 1000.0;
    if (elapsed > time_limit)
    {
      cout << "OOT\n";
      break;
    }
    cout << "\nLBs: ";
    for_each(d_LB_nodes.begin(), d_LB_nodes.end(), [](double val){ cout << val << ' '; });
    size_t node_idx = distance(d_LB_nodes.begin(), min_element(d_LB_nodes.begin(), d_LB_nodes.end()));
    cout << "\nExploring node " << node_idx << endl;
      // determine max_rounds based on tree_size
    bool branch = solve(types, node_idx, incumbent, tol, time_limit - elapsed, rcuts, fenchel, d_UB_global < 1e10 ? 0 : max_rounds);  // solve() also updates global bounds
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

  cout << "number of nodes: " << d_nodes.size() << '\n';
  cout << "tree size: " << tree_size << '\n';
  cout << "computation time: " << chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t1).count() / 1000.0 << '\n';
  return incumbent;
}

    
    






