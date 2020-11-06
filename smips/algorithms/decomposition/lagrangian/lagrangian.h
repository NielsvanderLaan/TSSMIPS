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
    vector<GRBVar> d_z_vars;
    vector<GRBVar> d_y_vars;
    GRBVar d_theta;
    bool d_rcut;

    Lagrangian(GRBEnv &env, Problem &problem, size_t s);
    Lagrangian(const Lagrangian &other);

    void update(vector<double> &pi);
    void add_cut(BendersCut &cut);
    void reverse_cut(double UB);
    void update_bound(size_t var, double val, bool lower);
    double solve();
    BendersCut strong_cut(vector<double> &pi);
    vector<double> z_vals();
};

#endif
