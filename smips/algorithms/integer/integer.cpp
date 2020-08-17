#include "integer.h"

bool is_integer(double val, double tol)
{
  return abs(val - floor(val + 0.5)) <= tol;
}