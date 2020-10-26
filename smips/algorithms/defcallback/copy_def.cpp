#include "DEF.h"

DEF::DEF(const DEF &other)
:
d_model(other.d_model)
{
  d_xvars = vector<GRBVar>(other.d_xvars.size());
  for (size_t var = 0; var != d_xvars.size(); ++var)
    d_xvars[var] = d_model.getVarByName("x_" + var);

  d_theta = d_model.getVarByName("theta");

}