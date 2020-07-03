#include "lagrangian.h"

Lagrangian::Lagrangian(const Lagrangian &other)
:
  d_model(other.d_model),
  d_n1(other.d_n1), 
  d_m2(other.d_m2), 
  d_n2(other.d_n2),
  d_problem(other.d_problem)
{
  d_constrs = d_model.getConstrs();      // heap allocated (deallocated in dtor)

  GRBVar *vars = d_model.getVars();      // heap allocated (deallocated at end of copy ctor)
  
  d_z_vars = new GRBVar[d_n1];           // heap allocated (deallocated in dtor)
  d_y_vars = new GRBVar[d_n2];

  copy_n(vars, d_n1, d_z_vars);
  copy_n(vars + d_n1, d_n2, d_y_vars);


  delete[] vars;                         // as announced
}