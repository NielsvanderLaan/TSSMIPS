#include "master.h"

double Master::violation()
{
  double violation, resid;
  GRBgetdblattr(d_cmodel, "ConstrVio", &violation);
  GRBgetdblattr(d_cmodel, "ConstrResidual", &resid);
  return violation + resid;
}