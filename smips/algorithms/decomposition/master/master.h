#ifndef MASTER_H
#define MASTER_H

#include "gurobi_c++.h"
#include "gurobi_c.h"
#include <vector>
#include <iostream>

#include "../../../problem_data/problem.h"
#include "../cut/benderscut.h"
#include "../../integer/integer.h"

using namespace std;

class Master
{
  public:  
    //GRBVar *d_xVars;
    //GRBVar d_theta;
    //GRBModel d_model;
    GRBmodel *d_cmodel;
      
      // internal storage: only valid for regular l-shaped and lbda cuts
      // slack variable identities s = kappa * theta + beta * x - gamma
    vector<double> d_kappa;
    vector<vector<double>> d_beta;
    vector<double> d_gamma;
    
    size_t d_n1;
    size_t d_p1;
    double d_L;
    size_t d_nSlacks;
    bool d_zk_safe;
    
    struct Solution
    {
      vector<double> xVals;
      double thetaVal;
      bool infeasible;
    };
    
    void update_bounds(size_t var, double val, bool lower); 
     
    Master(GRBEnv &env, GRBenv *c_env, Problem &problem, bool zk_safe); // initializes d_model and its variables
    Master(const Master &other);
    ~Master(); // deletes vars and frees d_cmodel
    
      // adds cut theta >= beta^T x + gamma, if this cut is violated (ret =  true), else cut is not added (ret = false). 
    bool add_ald_cut(double *beta, double gamma, double tau, double *x, double theta, double tol);
    bool add_cut(BendersCut cut, Solution sol, double tol);  // adds the cut kappa theta >= beta^T x + gamma 
    
    vector<BendersCut> round_of_cuts();
    BendersCut gmi_cut(size_t row, double a0);
      // auxiliary functions to derive gmi cut
    vector<double> extract_row(size_t row);
    void compute_cut(vector<double> &tab_row, double a0, double &coef_theta, double *coef_x);
    BendersCut transform_cut(double coef_theta, double *coef_x);
    
    Solution solve(); // solves the model 
};

#endif
