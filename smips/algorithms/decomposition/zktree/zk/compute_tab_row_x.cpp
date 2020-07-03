#include "zk.h"

void ZK::compute_tab_row_x(double *tab_row_x, int nVarsMaster, int row, GRBmodel *master)
{
  fill_n(tab_row_x, nVarsMaster, 0);
            
            // extracting i-th row of B^-1 (basis matrix of subproblem)
  int e_i_ind[1] = {row};  double e_i_val[1] = {1.0};
  GRBsvec e_i {1, e_i_ind, e_i_val};  // unit vector
  int Brow_ind[d_nConstrs];  double Brow_val[d_nConstrs];
  GRBsvec Brow {d_nConstrs, Brow_ind, Brow_val};   // result vector
  GRBBSolve(d_model, &e_i, &Brow);   // extracting ith row of B^-1

            // extracting master basis information
  int nConsMaster;                                      // number of constraints in master problem (including cuts)
  GRBgetintattr(master, "NumConstrs", &nConsMaster);
  int master_bhead[nConsMaster];
  GRBgetBasisHead(master, master_bhead);
  
            // computing (BW^-1)_i * T_{BA}
  double BinvTBA[nConsMaster];
  fill_n(BinvTBA, nConsMaster, 0);
  for (size_t idx = 0; idx != nConsMaster; ++ idx)
  {
    int basic_var = master_bhead[idx];  // column index of (tau, Tmat)    
    if (basic_var > d_n1)               // corresponds to slack variable, so 
      continue;                         // corresponding column of T is zero vector: skip 
          
    if  (basic_var == 0)            // then multiply by tau (column 0 of (tau, Tmat))
    {
      for (size_t nz = 0; nz != Brow.len; ++nz)     // loop over nonzeros of Brow
        BinvTBA[idx] += Brow.val[nz] * d_tau[Brow.ind[nz]];
    
    } else
    {
      for (size_t nz = 0; nz != Brow.len; ++nz)  // loop over nonzeros of Brow
        BinvTBA[idx] += Brow.val[nz] * d_Tmat[Brow.ind[nz]][basic_var - 1];    
    }     
  }

  for (size_t col = 0; col != nVarsMaster; ++col)
  {
    int BA_col_ind[nConsMaster];                            // extract column of master simplex tableau
    double BA_col_val[nConsMaster]; 
    GRBsvec BA_col{nConsMaster, BA_col_ind, BA_col_val};
    GRBBinvColj(master, col, &BA_col);


    for (size_t nz = 0; nz != BA_col.len; ++nz)
      tab_row_x[col] -= BinvTBA[BA_col.ind[nz]] * BA_col.val[nz];      
  }
  
      // add (B^W)^-1_i * tau to tab_row_x[0]
   for (size_t nz = 0; nz != Brow.len; ++nz)
    tab_row_x[0] += Brow.val[nz] * d_tau[Brow.ind[nz]];
      // add (B^W)^-1_i * Tmat to tab_row_x
  for (size_t idx = 0; idx != d_n1; ++idx) // no need to loop over slacks: corresponding rows of T are zero vectors
  {
    for (size_t nz = 0; nz != Brow.len; ++nz)
      tab_row_x[idx+1] += d_Tmat[Brow.ind[nz]][idx] * Brow.val[nz];
  }
}







