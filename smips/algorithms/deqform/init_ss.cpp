#include "deqform.h"

void DeqForm::init_ss(size_t n1, size_t n2, size_t p2, size_t m2, size_t S,
                               size_t ss_leq, size_t ss_geq,
                               double *lb, double *ub,
                               double *probs, vector<double> q,
                               vector<vector<double>> &Tmat, vector<vector<double>> &Wmat,
                               vector<vector<double>> &q_omega,
                               vector<vector<vector<double>>> &W_omega,
                               vector<vector<vector<double>>> &T_omega,
                               vector<vector<double>> &omega, bool fix_rec, bool fix_tech)
{
      // variable types    
  char vTypes2[n2];
  fill_n(vTypes2, p2, GRB_INTEGER);
  fill_n(vTypes2 + p2, n2 - p2, GRB_CONTINUOUS);

      //constraint senses
  char senses2[m2];
  fill(senses2,                   senses2 + ss_leq,          GRB_LESS_EQUAL);
  fill(senses2 + ss_leq,          senses2 + ss_leq + ss_geq, GRB_GREATER_EQUAL);
  fill(senses2 + ss_leq + ss_geq, senses2 + m2,              GRB_EQUAL);
  
      // for each scenario: add variables and constraints
  double prob;
  for (size_t s = 0; s != S; ++s)
  {
    vector<double> costs = fix_rec ? q : q_omega[s];
    for_each(costs.begin(), costs.end(), [probs, s](double &val){ val *= probs[s]; });

    double *rhsOmega = omega[s].data(); 
        // add variables
    GRBVar *yVars = d_model.addVars(lb, ub, costs.data(), vTypes2, NULL, n2);
    
        // lhs expression of second-stage constraints, Wy will be added in a loop
    vector<GRBLinExpr> TxWy(m2);
    vector<vector<double>> &tech = fix_tech ? Tmat : T_omega[s];
    vector<vector<double>> &rec_mat = fix_rec ? Wmat : W_omega[s];
    for (size_t conIdx = 0; conIdx != m2; ++conIdx)
    {
      TxWy[conIdx].addTerms(tech[conIdx].data(), d_xVars, n1);
      TxWy[conIdx].addTerms(rec_mat[conIdx].data(), yVars, n2);
    }


        // add constraints
    delete[] d_model.addConstrs(TxWy.data(), senses2, rhsOmega, NULL, m2);
    delete[] yVars;
  }   
}