#ifndef PSLP_H
#define PSLP_H

#include "../../../problem_data/problem.h"
#include "../master/master.h"
#include "../zktree/zk/zk.h"

class Pslp
{
  public:
    size_t d_n1, d_S;
    
    vector<ZK> d_zk;   
    vector<double> &d_probs;
    
    Pslp(GRBEnv &env, GRBenv *c_env, Problem &problem);
    Pslp(const Pslp &other);


    BendersCut best_zk_cut(Master::Solution sol, Master &master, bool lap_cuts = false, bool affine = false, size_t maxRounds = 25, double tol = 1e-4);
    
    void update_bounds(size_t var, double val, bool lower); // calls cglp::update_bound() via zk::update_bound()
    void add_cglp_rows(double *coef_x, double coef_theta, double *coef_y, double rhs);
    void reverse_cut(double UB);
};




#endif