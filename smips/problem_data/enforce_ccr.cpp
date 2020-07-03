#include "problem.h"

void Problem::enforce_ccr(double penalty)
{
  // we have to introduce new continuous variables (change n2)
  // we have to update q
  // we have to update W
  
  // we have to distuingish between <=, >=, and ==
  
  size_t ss_eq = d_m2 - d_ss_leq - d_ss_geq;
  size_t nArtVars = d_ss_leq + d_ss_geq + 2*ss_eq;   
  
  d_gen.append_r(d_q, penalty, nArtVars);      // updating costs
  d_gen.append_zeros(d_Wmat, nArtVars);        // appending columns (all-zeros)

  for (auto &q_vec : d_q_omega)
    d_gen.append_r(q_vec, penalty, nArtVars);
  for (auto &rec_mat : d_W_omega)
    d_gen.append_zeros(rec_mat, nArtVars);


  size_t col = d_n2;      
  size_t con = 0;        // constraint counter
  
  for (; con != d_ss_leq; ++con)
  {
    d_Wmat[con][col] = -1.0;
    for (auto &rec_mat : d_W_omega)
      rec_mat[con][col] = -1.0;
    ++col;
  }
  
  for (; con != d_ss_leq + d_ss_geq; ++con)
  {
    d_Wmat[con][col] = 1.0;
    for (auto &rec_mat : d_W_omega)
      rec_mat[con][col] = 1.0;
    ++col;
  }
  
  for (; con != d_m2; ++con)
  {
    d_Wmat[con][col] = -1.0;
    for (auto &rec_mat : d_W_omega)
      rec_mat[con][col] = -1.0;
    ++col;

    d_Wmat[con][col] = 1.0;
    for (auto &rec_mat : d_W_omega)
      rec_mat[con][col] = 1.0;
    ++col;
  }
  
  d_n2 += nArtVars;
  for (size_t var = 0; var != nArtVars; ++var)
  {
    d_l2.push_back(0);
    d_u2.push_back(GRB_INFINITY);
  }
  
}