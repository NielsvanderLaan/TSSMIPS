#ifndef LAGRANGIAN_H
#define LAGRANGIAN_H

#include "gurobi_c++.h"

#include "../../../problem_data/problem.h"

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

    void update(double *rhs, size_t s, double *pi);
    double solve();
};

#endif
