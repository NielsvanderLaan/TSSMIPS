#include "aggregator.h"

void Aggregator::init_vw(Problem &problem)
{
  vector<char> y_types(problem.d_n2, GRB_CONTINUOUS);
  fill_n(y_types.begin(), problem.d_p2, GRB_INTEGER);
  GRBVar *vw_vars = d_vw.addVars(problem.d_l2.data(), problem.d_u2.data(), problem.d_q.data(), y_types.data(), NULL, y_types.size());

  GRBLinExpr Wy[problem.d_m2];
  for (size_t con = 0; con != problem.d_m2; ++con)
    Wy[con].addTerms(problem.d_Wmat[con].data(), vw_vars, problem.d_n2);
  vector<double> rhs(problem.d_m2, 0.0);

  vector<char> senses(rhs.size(), GRB_EQUAL);
  fill_n(senses.begin(), problem.d_ss_leq, GRB_LESS_EQUAL);
  fill_n(senses.begin() + problem.d_ss_leq, problem.d_ss_geq, GRB_GREATER_EQUAL);

  delete[] d_vw.addConstrs(Wy, senses.data(), rhs.data(), NULL, rhs.size());
  delete[] vw_vars;

  d_vw.update();
}