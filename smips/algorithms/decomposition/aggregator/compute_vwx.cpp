#include "aggregator.h"

double Aggregator::compute_vwx(double *x, size_t s)
{
    // computing rhs
    
  vector<double> rhs(d_omega[s]);
  for (size_t con = 0; con != rhs.size(); ++con)
    rhs[con] -= inner_product(d_Tmat[con].begin(), d_Tmat[con].end(), x, 0.0); 
    // setting rhs and solving model
  GRBConstr *cons = d_vw.getConstrs();
  d_vw.set(GRB_DoubleAttr_RHS, cons, rhs.data(), rhs.size());
  d_vw.optimize();
    // clean up
  delete[] cons;
  return d_vw.get(GRB_DoubleAttr_ObjVal);  
}