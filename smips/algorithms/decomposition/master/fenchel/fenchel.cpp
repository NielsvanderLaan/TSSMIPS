#include "fenchel.h"

Fenchel::Fenchel(Problem &problem, GRBEnv &env)
:
  d_problem(problem),
  d_mp(env),
  d_sub(env),
  d_rcut(false)
{
  // Initializing mp
  double M = 1.0;
  vector<double> lb(problem.d_n1, -M);
  vector<double> ub(problem.d_n1, M);

  d_alpha = d_mp.addVar(-M, M, -1.0, GRB_CONTINUOUS, "alpha");
  GRBVar *beta = d_mp.addVars(lb.data(), ub.data(), NULL, NULL, NULL, problem.d_n1);
  d_beta = vector<GRBVar> (beta, beta + problem.d_n1);
  delete[] beta;
  d_kappa = d_mp.addVar(-M, M, 0.0, GRB_CONTINUOUS, "kappa");

  // Initializing sub
  vector<char> vtypes(problem.d_n1, GRB_CONTINUOUS);
  fill_n(vtypes.begin(), problem.d_p1, GRB_INTEGER);
  GRBVar *xvars = d_sub.addVars(problem.d_l1.data(), problem.d_u1.data(), NULL, vtypes.data(), NULL, problem.d_n1);
  d_xvars = vector<GRBVar>(xvars, xvars + problem.d_n1);
  delete[] xvars;

  d_theta = d_sub.addVar(problem.d_L, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "theta");

  // imposing Ax ~ b
  vector<char> senses(problem.d_m1, GRB_EQUAL);
  fill_n(senses.begin(), problem.d_fs_leq, GRB_LESS_EQUAL);
  fill_n(senses.begin() + problem.d_fs_leq, problem.d_fs_geq, GRB_GREATER_EQUAL);
  GRBLinExpr Ax[problem.d_m1];
  for (size_t con = 0; con != problem.d_m1; ++con)
    Ax[con].addTerms(problem.d_Amat[con].data(), d_xvars.data(), d_xvars.size());
  delete[] d_sub.addConstrs(Ax, senses.data(), problem.d_b.data(), NULL, senses.size());

  d_mp.update();
  d_sub.update();
}