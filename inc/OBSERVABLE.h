// ----------------------------------------------------------------------------- 
// Author:  James Lambert
// Created: February 2nd, 2018
//
// Template class for SSE observables. Useful for keeping track of errors. 
// -----------------------------------------------------------------------------

#ifndef OBSERVABLE_H
#define OBSERVABLE_H

#include <vector>

namespace SSE
{
  template <class type>
  class observable
  {
    private:

      int               _bs;      // bin size for the observable
      int               _ind;     // current index in bin size
      double            _ave;     // average for observable
      double            _err;     // error in average

      type              _sum;     // total value for current bin
      type              _sq;      // total sq for current bin
      
      std::vector<int>  _binave;
      std::vector<int>  _binerr;

    public:

      observable(const int& b) : _bs(b), _ind(0), _ave(0.0), _err(0.0){}
      
      void push_back(type);

      double ave();
      double err();
  };

  void observable::push_back(const type& val)
  {
    if(index = 0)
    {
      _sum = val;
      _sq  = val * val;
      index += 1;
    }
    else if(index < _bs-1)
    {
      _sum += val;
      _sq  += val * val;
      index += 1;
    }
    else
    { 
      double a = (double)_sum / (double)_bs;
      double e = a*a - ((double)sq/(double)_bs);
      
      _binave.push_back(a);
      _binerr.push_back(e);
      
      index = 0;
    }
  }
  
  double ave() const
  {
    double s = 0.0;
    for(auto it=_binave.begin(); it!=_binave.end(); it++)
      s += *it;
    return s / (double)_binave.size();
  }

  double err() const
  {
    double s = 0.0;
    for(auto it=_binerr.begin(); it!=binave.end(); it++)
      s += *it;
    return s / (double)_binerr.size();
  }

}

#endif // OBSERVABLE_H


