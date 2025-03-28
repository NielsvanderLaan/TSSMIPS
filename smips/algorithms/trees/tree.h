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
    
    Tree(GRBEnv &env, GRBenv *c_env, Problem &problem, vector<Type> &types);
    Tree(GRBEnv &emv, GRBenv *c_env, Problem &problem, Benders &root);
    ~Tree();
    
    vector<double> bab(vector<Type> types, bool rcuts = true, bool fenchel = true, int max_rounds = 5, double tol = 1e-4, double time_limit = 24*3600);
        // auxiliary functions
    bool solve(vector<Type> types, size_t node_idx, vector<double> &incumbent, double local_tol, double time_limit, bool rcuts, bool fenchel, int max_rounds);
    
    struct Split
    {
      int var;
      double left, right;  
    };
    bool add_branch(size_t node_idx, Split split);
    Split branch_var(size_t node_idx);
    void fathom();

    int c_branch_var(Benders *node, double *x);
    int c_branch_var_diam(Benders* node);
    int i_branch_var(double *x);
};



#endif