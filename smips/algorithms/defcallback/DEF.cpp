#include "DEF.h"

DEF::DEF(Problem &problem, GRBEnv &env)
:
  d_model(env)
{
  d_model.set(GRB_IntParam_OutputFlag, 1);
  d_model.set(GRB_IntParam_LazyConstraints, 1);

  char vTypes[problem.d_n1];
  fill_n(vTypes, problem.d_p1, GRB_INTEGER);
  fill_n(vTypes + problem.d_p1, problem.d_n1 - problem.d_p1, GRB_CONTINUOUS);

  vector<string> names(problem.d_n1);
  for (size_t var = 0; var != problem.d_n1; ++var)
    names[var] = "x_" + var;
  GRBVar *xvars = d_model.addVars(problem.d_l1.data(), problem.d_u1.data(), problem.d_c.data(), vTypes, names.data(), problem.d_n1);
  d_xvars = vector<GRBVar>(xvars, xvars + problem.d_n1);
  delete[] xvars;
  //vector<int> ones(d_xvars.size(), 1);
  //d_model.set(GRB_IntAttr_BranchPriority, d_xvars.data(), ones.data(), ones.size());


  GRBLinExpr Ax[problem.d_m1];
  for (size_t con = 0; con!= problem.d_m1; ++con)
    Ax[con].addTerms(problem.d_Amat[con].data(), d_xvars.data(), d_xvars.size());
  vector<char> senses(problem.d_m1, GRB_EQUAL);
  fill_n(senses.begin(),                    problem.d_fs_leq, GRB_LESS_EQUAL);
  fill_n(senses.begin() + problem.d_fs_leq, problem.d_fs_geq, GRB_GREATER_EQUAL);

  delete[] d_model.addConstrs(Ax, senses.data(), problem.d_b.data(), NULL, senses.size());

  d_theta = d_model.addVar(problem.d_L, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "theta");

  char ytypes[problem.d_n2];
  fill_n(ytypes, problem.d_p2, GRB_INTEGER);
  fill_n(ytypes + problem.d_p2, problem.d_n2 - problem.d_p2, GRB_CONTINUOUS);

  vector<char> senses2(problem.d_m2, GRB_EQUAL);
  fill_n(senses2.begin(),                    problem.d_ss_leq, GRB_LESS_EQUAL);
  fill_n(senses2.begin() + problem.d_ss_leq, problem.d_ss_geq, GRB_GREATER_EQUAL);

  GRBLinExpr qy;
  for (size_t s = 0; s != problem.d_S; ++s)
  {
    GRBVar *yvars = d_model.addVars(problem.d_l2.data(), problem.d_u2.data(), NULL, ytypes, NULL, problem.d_n2);
          // adding Wy + Tx = h
    vector<vector<double>> &rm = problem.d_fix_rec ? problem.d_Wmat : problem.d_W_omega[s];
    GRBLinExpr WyTx[problem.d_m2];
    for (size_t con = 0; con != problem.d_m2; ++con)
    {
      WyTx[con].addTerms(rm[con].data(), yvars, rm[con].size());
      WyTx[con].addTerms(problem.d_Tmat[con].data(), d_xvars.data(), d_xvars.size());
    }
    delete[] d_model.addConstrs(WyTx, senses2.data(), problem.d_omega[s].data(), NULL, senses2.size());

    vector<double> q = problem.d_fix_rec ? problem.d_q : problem.d_q_omega[s];
    double prob = problem.d_probs[s];
    for_each(q.begin(), q.end(), [prob](double &val){ val *= prob; });

    qy.addTerms(q.data(), yvars, q.size());
    delete[] yvars;
  }
  d_model.addConstr(d_theta == qy);
  d_model.update();
}