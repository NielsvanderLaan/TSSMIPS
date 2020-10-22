#ifndef TSSMIPS_DEF_H
#define TSSMIPS_DEF_H


#include <vector>
#include "gurobi_c++.h"
#include "../../problem_data/problem.h"

class DEF
{
  public:
    DEF(Problem &problem, GRBEnv &env);
    DEF(const DEF &other) = delete;
    GRBModel d_model;
    vector<GRBVar> d_xvars;
    GRBVar d_theta;

    void solve();
};

#endif //TSSMIPS_DEF_H
