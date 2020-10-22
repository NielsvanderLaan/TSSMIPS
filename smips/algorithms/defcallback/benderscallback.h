#ifndef TSSMIPS_BENDERSCALLBACK_H
#define TSSMIPS_BENDERSCALLBACK_H

#include <vector>
#include "gurobi_c++.h"
#include "gurobi_c.h"
#include "../../problem_data/problem.h"
#include "../decomposition/benders.h"
#include "DEF.h"

using namespace std;

class BendersCallback: public GRBCallback
{
  public:
    vector<GRBVar> d_xvars;
    GRBVar d_theta;
    Benders d_ben;
    Problem &d_problem;
    int d_ncuts;

    BendersCallback(Problem &problem, GRBEnv &env, GRBenv* c_env, DEF &def);
    BendersCallback(const BendersCallback &other) = delete;

  protected:
    void callback();

};

#endif //TSSMIPS_BENDERSCALLBACK_H
