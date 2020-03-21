#ifndef ZKTREE_H
#define ZKTREE_H

#include "zk/zk.h"


class ZkTree
{
  public:
    vector<ZK *> d_nodes;  // LP-relaxations of nodes

    
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
    
    // generate cut
    // add_row_to_cglp_subproplems();    // all or one?
     
};

#endif