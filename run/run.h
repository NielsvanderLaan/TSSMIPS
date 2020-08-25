#ifndef TSSMIPS_RUN_H
#define TSSMIPS_RUN_H

#include <iostream>
#include <chrono>
#include <vector>

#include "gurobi_c++.h"
#include "gurobi_c.h"

#include "problem.h"
#include "deqform.h"
#include "benders.h"
#include "tree.h"

void run_ssv_ld_gaps(Data &rand, GRBEnv &env, GRBenv *c_env);
void lbda_scheme(Problem &problem, GRBEnv &env, GRBenv *c_env, size_t nIter);
void solve_dcap(Data &rand, GRBEnv &env, GRBenv *c_env);

#endif //TSSMIPS_RUN_H
