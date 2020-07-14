#include "benders.h"

Benders::~Benders()
{
  delete[] d_xvals;
  delete[] d_incumbent;
}