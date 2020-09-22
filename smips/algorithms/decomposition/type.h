#ifndef TSSMIPS_TYPE_H
#define TSSMIPS_TYPE_H

#include <string>
#include <vector>

enum Type {LP, SB, LR, SC_CPT, SC_BAB, SC_RG, LBDA};

static std::string name(Type type)
{
  switch (type)
  {
    case LP:      return "LP cut";
    case SB:      return "SB cut";
    case LR:      return "LR cut";
    case SC_CPT:  return "CPT scaled cut";
    case SC_BAB:  return "BAB scaled cut";
    case SC_RG:   return "RG scaled cut";
    case LBDA:    return "LBDA cut";
    default:      return "unknown cut type";
  }
}

static std::vector<Type> string_to_type(char *types[], int nTypes)
{
  std::vector<Type> ret;
  for (size_t idx = 0; idx != nTypes; ++idx)
  {
    std::string type_string(types[idx]);
    if (type_string == "LP")      ret.push_back(LP);
    if (type_string == "SB")      ret.push_back(SB);
    if (type_string == "LR")      ret.push_back(LR);
    if (type_string == "ZK")      ret.push_back(SC_CPT);
    if (type_string == "CPT")     ret.push_back(SC_CPT);
    if (type_string == "BAB")     ret.push_back(SC_BAB);
    if (type_string == "RG")      ret.push_back(SC_RG);
    if (type_string == "LBDA")    ret.push_back(LBDA);
  }
  return ret;
}



#endif //TSSMIPS_TYPE_H
