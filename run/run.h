#ifndef TSSMIPS_RUN_H
#define TSSMIPS_RUN_H

#include <iostream>
#include <chrono>
#include <vector>

#include "gurobi_c++.h"
#include "gurobi_c.h"

#include "../smips/problem_data/problem.h"
#include "../smips/algorithms/deqform/deqform.h"
#include "../smips/algorithms/decomposition/benders.h"
#include "../smips/algorithms/trees/tree.h"

void run_ssv_ld_gaps(Data &rand, GRBEnv &env, GRBenv *c_env);
void lbda_scheme(Problem &problem, GRBEnv &env, GRBenv *c_env, size_t nIter);
void solve_dcap(bool lp_cuts, bool sb_cuts, bool zk_cuts, bool strong_cuts, Data &rand, GRBEnv &env, GRBenv *c_env);
void solve_sizes(bool lp_cuts, bool sb_cuts, bool zk_cuts, bool strong_cuts, Data &rand, GRBEnv &env, GRBenv *c_env);

#endif //TSSMIPS_RUN_H
