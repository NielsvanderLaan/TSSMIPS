#ifndef TREE_H
#define TREE_H

#include<algorithm>
#include "../decomposition/benders.h"

using namespace std;

class Tree
{
  public:
    Problem &d_problem;
    GRBenv *d_c_env;
    GRBEnv &d_env;
    
    vector<Benders*> d_nodes;
    
    double d_UB_global;
    double d_LB_global;  
    vector<double> d_LB_nodes;
    
    Tree(GRBEnv &env, GRBenv *c_env, Problem &problem);
    Tree(GRBEnv &emv, GRBenv *c_env, Problem &problem, Benders &root);
    ~Tree();
    
    vector<double> bab(bool affine = false, double tol = 1e-4);
        // auxiliary functions
    bool solve(size_t node_idx, vector<double> &incumbent, bool affine, double global_tol, double local_tol, int max_iter);
    
    struct Split
    {
      int var;
      double left, right;  
    };
    bool branch(size_t node_idx, Split split);
    Split branch_var(size_t node_idx);
    int c_branch_var(Benders *node, double *x);
    int i_branch_var(double *x);
};



#endif