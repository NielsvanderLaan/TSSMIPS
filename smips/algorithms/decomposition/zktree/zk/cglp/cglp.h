#ifndef CGLP_H
#define CGLP_H

#include <algorithm>
#include <iostream>

#include "gurobi_c.h"

#include "../../../../../problem_data/problem.h"
#include "../cut/cut.h"

class Cglp
{
  public:
    GRBModel d_model;
    
    size_t d_n1, d_p1, d_m1, d_n2, d_p2, d_m2, d_nMults;
    
    vector<GRBConstr> d_constrs1;    // corresponding to first term of disjunction
    vector<GRBConstr> d_constrs2;    // corresponding to second term of disjunction
    
    vector<int> d_l1_mults, d_u1_mults, d_l2_mults, d_u2_mults;
    
    // cuts are of the form Wrow^T y + Trow^T x + r theta >= h
    vector<GRBVar> d_Trow, d_Wrow;   // lengths: n1 and n2
    vector<GRBVar> d_lambda1, d_lambda2;  // one lambda for each term of the disjunction, initial length: m1 + m2 + 1
    GRBVar d_r, d_h;
    
    Cglp(Problem &problem, GRBEnv &env, size_t scenario);
    Cglp(const Cglp &other);
    Cglp(Cglp &&other) = delete;
    
    void create_disjunction(size_t var_idx, double val);  // creates the disjunction: y[k] <= val, y[k] >= val + 1
    
    void add_row(double *coef_x, double coef_theta, double *coef_y, double rhs);
    void update_bound(size_t var, double val, bool lower, bool fs);    // updates cglp, called if a  bound changes
    void set_obj(double *x, double theta, double *y);

    Cut generate_cut(double *x, double theta, double *y, size_t var_idx, double val);
};

#endif