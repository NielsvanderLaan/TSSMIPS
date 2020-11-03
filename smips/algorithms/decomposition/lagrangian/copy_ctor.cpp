#include "lagrangian.h"


Lagrangian::Lagrangian(const Lagrangian &other)
:
  d_model(nullptr),
  d_n1(other.d_n1), 
  d_m2(other.d_m2), 
  d_n2(other.d_n2),
  d_problem(other.d_problem),
  d_z_vars(other.d_n1),
  d_y_vars(other.d_n2),
  d_rcut(other.d_rcut)
{
/*
#pragma omp critical
{


  d_model.write("copy.lp");
  other.d_model.write("original.lp");

  exit(1);
}
*/

#pragma omp critical
  d_model = new GRBModel(*other.d_model);

  GRBConstr *cons = d_model->getConstrs();      // heap allocated (deallocated in dtor)
  d_constrs = vector<GRBConstr>(cons + d_problem.d_m1, cons + d_problem.d_m1 + d_problem.d_m2);
  delete[] cons;

  GRBVar *vars = d_model->getVars();      // heap allocated (deallocated at end of copy ctor)

  copy_n(vars, d_n1, d_z_vars.begin());
  copy_n(vars + d_n1, d_n2, d_y_vars.begin());
  d_theta = vars[d_n1 + d_n2];


  delete[] vars;                         // as announced

}