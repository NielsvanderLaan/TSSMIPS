#include "master.h"

// write round_of_cuts(), which returns row indices and a0

BendersCut Master::gmi_cut(size_t row, double a0)
{
        // extracting tableau row
  vector<double> tab_row = extract_row(row);
    
       // computing gomory cut coefficients
  double coef_theta;
  double coef_x[tab_row.size() - 1];               // original x variables and slacks     

  compute_cut(tab_row, a0, coef_theta, coef_x);    // return by argument
  
  return transform_cut(coef_theta, coef_x);
}