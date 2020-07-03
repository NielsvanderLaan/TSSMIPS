#include "problem.h"

double Problem::evaluate(double *x)
{
    // computing cx
  double cx = 0.0;
  for (size_t var = 0; var != d_n1; ++var)
    cx += d_c[var] * x[var];
    
    // computing z = Tx
  double Tx[d_m2];  
  for (size_t zvar = 0; zvar!= d_m2; ++zvar)
  {
    Tx[zvar] = 0.0; 
    for (size_t xvar = 0; xvar != d_n1; ++xvar)
      Tx[zvar] += d_Tmat[zvar][xvar] * x[xvar];
  }

  double Q = 0.0;

  GRBModel *sub = init_sub();
  for (size_t s = 0; s != d_S; ++s)
  {
    double *ws = d_omega[s].data();
    double rhs[d_m2];
    for (size_t zvar = 0; zvar != d_m2; ++zvar)
      rhs[zvar] = ws[zvar] - Tx[zvar];
    GRBConstr *constrs = sub->getConstrs();
    sub->set(GRB_DoubleAttr_RHS, constrs, rhs, d_m2);

    if (not d_fix_rec)
    {
      GRBVar *vars = sub->getVars();
      sub->set(GRB_DoubleAttr_Obj, vars, d_q_omega[s].data(), d_q_omega[s].size());
      vector<vector<double>> &Ww = d_W_omega[s];
      for (size_t row = 0; row != Ww.size(); ++row)
      {
        vector<GRBConstr> con_vec(Ww[row].size(), constrs[row]);
        sub->chgCoeffs(con_vec.data(), vars, Ww[row].data(), Ww[row].size());
      }
      delete[] vars;
    }

    delete[] constrs;

    sub->optimize();
    Q += d_probs[s] * sub->get(GRB_DoubleAttr_ObjVal);
  }
  delete sub;
  return cx + Q;   
}