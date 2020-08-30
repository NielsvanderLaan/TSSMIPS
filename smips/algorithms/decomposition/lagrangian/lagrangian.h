#ifndef LAGRANGIAN_H
#define LAGRANGIAN_H

#include "gurobi_c++.h"

#include "../../../problem_data/problem.h"
#include "../cut/benderscut.h"

using namespace std;

class Lagrangian
{
  public:
    GRBModel d_model;
    size_t d_n1, d_m2, d_n2;
    Problem &d_problem;
    vector<GRBConstr> d_constrs;
    GRBVar *d_z_vars;
    GRBVar *d_y_vars;
    Lagrangian(GRBEnv &env, Problem &problem);
    Lagrangian(const Lagrangian &other);
    ~Lagrangian();

    void update(size_t s,  vector<double> &pi);
    void update_pi( vector<double> &pi);
    double solve();
    BendersCut strong_cut(size_t s, vector<double> &pi);
    BendersCut lr_cut(size_t s, double *x, double vwx,  vector<double> &pi);
    vector<double> z_vals();
};

#endif
