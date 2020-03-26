#ifndef BendersCut_H
#define BendersCut_H

#include <algorithm>

using namespace std;

struct BendersCut
{
  double d_alpha;
  vector<double> d_beta;
  double d_tau;
  
  BendersCut operator*(double scale)
  {
    vector<double> beta(this->d_beta);
    for_each(beta.begin(), beta.end(), [scale](double &val){ val *= scale; });
    return BendersCut{ this->d_alpha * scale, beta, this->d_tau * scale };
  }
  
  BendersCut& operator+=(const BendersCut &right)
  {
    this->d_alpha += right.d_alpha;
    this->d_tau += right.d_tau;
    transform(this->d_beta.begin(), this->d_beta.end(), right.d_beta.begin(), this->d_beta.begin(), plus<double>());
    
    return *this;
  }
  
  
};

#endif