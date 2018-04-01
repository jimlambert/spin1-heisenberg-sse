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
      double            _err;     // error in average

      type              _sum;     // total value for current bin
      type              _sq;      // total sq for current bin
      
      std::vector<double>  _binave;
      std::vector<double>  _binerr;

    public:

      OBSERVABLE() : _bs(1000), _ind(0), _ave(0.0), _err(0.0){}
      
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
      _sq  = val * val;
      _ind += 1;
    }
    else if(_ind < _bs-1)
    {
      _sum += val;
      _sq  += val * val;
      _ind += 1;
    }
    else
    {
      _sum += val;
      _sq  += val*val; 
      double a = (double)_sum / (double)_bs;
      double e = std::sqrt((((double)_sq/(double)_bs) - a*a)/(double)_bs);
      _binave.push_back(a);
      _binerr.push_back(e);
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
    for(auto it=_binerr.begin(); it!=_binerr.end(); it++)
      s += *it;
 
    return s / (double)_binerr.size();
  }

}

#endif // OBSERVABLE_H


