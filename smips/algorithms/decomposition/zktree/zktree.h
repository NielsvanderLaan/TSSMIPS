#ifndef ZKTREE_H
#define ZKTREE_H

#include "zk/zk.h"


class ZkTree
{
  public:
    Problem &d_problem;
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
    vector<int> d_rcut_inds;
    double d_cglp_val;
    BendersCut d_candidate;
    
    ZkTree(GRBenv *env, GRBEnv &cpp_env, Problem &problem, size_t scenario);
    ZkTree(const ZkTree &other);     // copy ctor

    ~ZkTree();      // update dtor to delete d_subs objects
    
    void branch_cut(double *x, double theta, double rho, Master &master, bool cuts, size_t maxRounds, double tol);
    
        // auxiliary functions
              // solve node and updates global bounds accordingly
    void update_global_bounds(size_t node_idx, vector<double> &lb_nodes, double &UB, double *x, double theta, double rho, Master &master, bool cuts, size_t maxRounds, double tol);
              
              // calls node->solve() or node->optimize() and updates cglp, returns true if node is feasible 
    bool solve(size_t node_idx, double *x, double theta, double rho, Master &master, bool cuts, size_t maxRounds, double tol);
    void branch(size_t node_idx, vector<double> &lb_nodes, double &UB, double *x, double theta, double rho, Master &master, bool cuts, size_t maxRounds, double tol);
    void add_child(size_t node_idx);            // copies d_nodes[node_idx] and appends it to d_nodes and updates cglp
    void update_bound(size_t node_idx, size_t var_idx, bool lower);  // imposes branching restrictions
    bool is_feasible(ZK *node);                               // checks if lp solution of node satisfies integer requirements 
    size_t branch_var_idx(size_t node_idx, bool strong = true);      // branching rule   

    void add_row_to_cglp(const double *coeff_x, double coeff_theta, double coeff_eta, double rhs, size_t node_idx);
    void add_benders_cut(const BendersCut &cut); // calls add_row_to_cglp for every node_idx
    void update_fs_bounds(size_t var, double val, bool lower); // updates cglp to incorporate changes to fs bounds (resulting from branching)
    void reverse_cut(double UB);

    BendersCut generate_cut(double *x, double theta, double rho, Master &master, bool cuts, size_t maxRounds, double tol); // calls branch_cut() and computes cut by solving cglp

    void solve_cglp(double M = 1e8);
    void optimize();
    double max_coeff(BendersCut &cut);
    void cglp_bounds(double M, bool set);

    double cglp_val();
};

#endif