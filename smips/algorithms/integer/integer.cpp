#include "integer.h"

bool is_integer(double val)
{
  return abs(val - floor(val + 0.5)) <= 1e-5;
}