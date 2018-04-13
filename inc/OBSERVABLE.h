// ----------------------------------------------------------------------------- 
// Author:  James Lambert
// Created: February 2nd, 2018
//
// Template class for SSE observables. Useful for keeping track of errors. 
// -----------------------------------------------------------------------------

#ifndef OBSERVABLE_H
#define OBSERVABLE_H

#include <iostream>
#include <vector>
#include <cmath>

namespace SSE
{
  template <class type>
  class OBSERVABLE
  {
    private:

      int               _bs;      // bin size for the OBSERVABLE
      int               _ind;     // current _ind in bin size
      double            _ave;     // average for OBSERVABLE

      type              _sum;     // total value for current bin
      
      std::vector<double>  _binave;

    public:

      OBSERVABLE() : _bs(1000), _ind(0), _ave(0.0){}
      
      void push_back(const type&);

      double ave();
      double err();
  };

  template <class type>
  void OBSERVABLE<type>::push_back(const type& val)
  {
    if(_ind == 0)
    {
      _sum = val;
      _ind += 1;
    }
    else if(_ind < _bs-1)
    {
      _sum += val;
      _ind += 1;
    }
    else
    {
      _sum += val;
      double a = (double)_sum / (double)_bs;
      _binave.push_back(a);
      _ind = 0;
    }
  }
  
  template <class type>
  double OBSERVABLE<type>::ave()
  {
    double s = 0.0;
    for(auto it=_binave.begin(); it!=_binave.end(); it++)
      s += *it; 
    return s / (double)_binave.size();
  }

  template <class type>
  double OBSERVABLE<type>::err()
  {
    double s = 0.0;
    for(auto it=_binave.begin(); it!=_binave.end(); it++)
      s += *it;
    s = s / (double)_binave.size();
    double tot = 0.0;
    for(auto it=_binave.begin(); it!=_binave.end(); it++)
      tot += (*it - s)*(*it-s);
    int B = _binave.size(); 
    return std::sqrt(tot / ((double)B*((double)B-1)));
  }

}

#endif // OBSERVABLE_H


