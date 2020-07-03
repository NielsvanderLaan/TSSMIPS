#include "problem.h"

GRBModel* Problem::init_sub()
{
  GRBModel *sub = new GRBModel(d_env);
  char vTypes[d_n2];
  fill(vTypes, vTypes + d_p2, GRB_INTEGER);  
  fill(vTypes + d_p2, vTypes + d_n2, GRB_CONTINUOUS);
      // add variables

  GRBVar *vars = sub->addVars(d_l2.data(), d_u2.data(), d_q.data(), vTypes, NULL, d_n2);
  
      // constraint senses
  char senses[d_m2];
  fill(senses,                       senses + d_ss_leq,            GRB_LESS_EQUAL);
  fill(senses + d_ss_leq,            senses + d_ss_leq + d_ss_geq, GRB_GREATER_EQUAL);
  fill(senses + d_ss_leq + d_ss_geq, senses + d_m2,                GRB_EQUAL);

  GRBLinExpr Wy[d_m2];
  for (size_t conIdx = 0; conIdx != d_m2; ++conIdx)
    Wy[conIdx].addTerms(d_Wmat[conIdx].data(), vars, d_n2);

      // add constraints
  vector<double> rhs(d_m2);
  delete[] sub->addConstrs(Wy, senses, rhs.data(), NULL, d_m2);
  delete[] vars;
  sub->update();
  return sub;
}