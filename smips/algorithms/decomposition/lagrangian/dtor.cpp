#include "lagrangian.h"

Lagrangian::~Lagrangian()
{
  if (d_model != nullptr)
    delete d_model;
}