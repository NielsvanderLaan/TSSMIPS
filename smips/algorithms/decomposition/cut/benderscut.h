#ifndef BendersCut_H
#define BendersCut_H

#include <algorithm>

using namespace std;

struct BendersCut
{
  double d_alpha;
  vector<double> d_beta;
  double d_tau;

  bool d_feas_cut;

  BendersCut()
  :
    d_feas_cut(false)
  {}

  BendersCut(double alpha, vector<double> beta, double tau, bool feas_cut = false)
  :
    d_alpha(alpha),
    d_beta(beta),
    d_tau(tau),
    d_feas_cut(feas_cut)
  {}

  
  BendersCut operator*(double scale)
  {
    vector<double> beta(this->d_beta);
    for_each(beta.begin(), beta.end(), [scale](double &val){ val *= scale; });
    return BendersCut{ this->d_alpha * scale, beta, this->d_tau * scale, this->d_feas_cut};
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