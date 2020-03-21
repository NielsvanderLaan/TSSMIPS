#ifndef ZK_H
#define ZK_H

#include <algorithm>
#include <iostream>
#include <string>

#include "gurobi_c++.h"
#include "gurobi_c.h"

#include <math.h>

#include "../../../../problem_data/problem.h"
#include "../../master/master.h"
#include "../../../integer/integer.h"
#include "cglp/cglp.h"
#include "cut/cut.h"
#include "../../cut/benderscut.h"


using namespace std;

class ZK
{
  public:
    size_t d_n1, d_p1, d_n2, d_p2, d_m2;
    int d_nVars, d_nConstrs;
    
    Cglp d_cglp;
    
    // constraints are of the form Wy ~ omega - T * x - tau * theta
    // includes constraints, upper and lower bounds, and cutting planes
    // x and y do not include slacks (constraints do not feature slacks).
    // we keep track for computing the rhs, for subgradient inequalities, and for transforming the cuts (so that they do not feature slacks)
    vector<vector<double>> d_Wmat;
    vector<vector<double>> d_Tmat;
    vector<double> d_tau;
    vector<double> d_omega;
    vector<int> d_signs;              // type of inequality -1 means <= , 0 means ==, and 1 means >=
    vector<int> d_lb_inds, d_ub_inds; // constraint indices of upper and lower bounds (-1 if no bounds are imposed (lb = 0, or ub = 1e20))
    
    
    double d_objVal;
    vector<double> d_yvals;
    GRBmodel *d_model;
    
    ZK(GRBenv *env, GRBEnv &cpp_env, Problem &problem, size_t scenario);
    ZK(const ZK &other);
    ZK(ZK &&other);
    ~ZK();
    
  
    void update(double *x, double theta);    // computes and updates rhs
    bool solve(double *x, double theta, Master &master, size_t maxRounds, bool gomory = true, double tol = 1e-2);    // returns false if model is infeasible (may happen due to branching)
    BendersCut subgradient();  // v_w(x) >= alpha + beta^T x + tau * theta  
    bool optimize();    // returns false if model is infeasible (may happen due to branching)
    
        // computes the rhs (coef_rhs - coef_theta * theta - coef_x^T x), adds the cut,
        // stores cut coeffients (omega, Tmat, and tau) for future calls to update(), and calls add_cglp_row()
    bool add_cut(Cut cut, double *x, double theta, double tol, size_t conIdx);                                  
    void update_bound(int var, double val, bool lower, bool fs);                     // updates lower/upper bounds in first/second-stage
    void add_cglp_row(double *coef_x, double coef_theta, double *coef_y, double rhs);   // calls cglp::add_row() to add cutting plane or optimality cut to cglp
  
    Cut generate_gmi_cut(Master &master, size_t row, double yval);    // generates a gmi cut
    void compute_tab_row_x(double *tab_row_x, int nVarsMaster, int row, GRBmodel *master);    // support functions for generate_gmi_cut()
    void compute_tab_row_y(double *tab_row_y, int row);                                       // idem
    void gmi_cut(double *tab_row_x, double *tab_row_y, double a0, double *coef_x, double *coef_y, double &coef_theta, int nVarsMaster);  // idem
    void transform_cut(double *coef_x, double *coef_y, double &coef_theta, double &coef_rhs, vector<double> &kappa, vector<vector<double>> &beta, vector<double> &gamma, size_t nSlacks);  // idem  
    
    double probe(size_t var_idx, double val, bool lower); // temporarily adjusts lb/ub of y[var_idx] to val, and returns increase in objective value. On return, leaves object unchanged.
  
 
};


#endif









