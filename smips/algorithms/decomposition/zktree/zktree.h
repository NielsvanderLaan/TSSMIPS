#ifndef ZKTREE_H
#define ZKTREE_H

#include "zk/zk.h"


class ZkTree
{
  public:
    vector<ZK *> d_nodes;  // LP-relaxations of nodes
    double d_L;
    GRBModel d_cglp;
        // kappa eta + beta^T x + tau theta >= alpha
    GRBVar d_alpha;
    vector<GRBVar> d_beta;
    GRBVar d_tau;
    GRBVar d_kappa;
    vector<vector<GRBConstr>> d_constrs; // for each term in the disjunction, we store the constraints (n + 3)
    vector<vector<GRBVar>> d_lambda;
    vector<vector<int>> d_lb_mult_inds, d_ub_mult_inds;
    
    ZkTree(GRBenv *env, GRBEnv &cpp_env, Problem &problem, size_t scenario);
    // copy ctor
    ~ZkTree();      // update dtor to delete d_subs objects
    
    void branch_cut(double *x, double theta, Master &master, size_t maxRounds, bool gomory = true);  
    
        // auxiliary functions
              // solve node and updates global bounds accordingly
    void update_global_bounds(size_t node_idx, vector<double> &lb_nodes, double &UB, double *x, double theta, Master &master, size_t maxRounds, bool gomory); 
              
              // calls node->solve() or node->optimize() and updates cglp, returns true if node is feasible 
    bool solve(size_t node_idx, double *x, double theta, Master &master, size_t maxRounds, bool gomory);     
    void branch(size_t node_idx, vector<double> &lb_nodes, double &UB, double *x, double theta, Master &master, size_t maxRounds, bool gomory);     
    void add_child(size_t node_idx);            // copies d_nodes[node_idx] and appends it to d_nodes and updates cglp
    void update_bound(size_t node_idx, size_t var_idx, bool lower);  // imposes branching restrictions and updates cglp
    bool is_feasible(ZK *node);                               // checks if lp solution of node satisfies integer requirements 
    size_t branch_var_idx(size_t node_idx, bool strong = true);      // branching rule   

    void add_row_to_cglp(double *coeff_x, double coeff_theta, double coeff_eta, double rhs, size_t node_idx);
    void add_benders_cut(double *coeff_x, double coeff_theta, double rhs); // calls add_row_to_cglp for every node_idx
    void update_fs_bounds(size_t var, double val, bool lower); // updates cglp to incorporate changes to fs bounds (resulting from branching)

    BendersCut generate_cut(double *x, double theta); // IMPORTANT: take into account d_L

    // generate cut
     
};

#endif