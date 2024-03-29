// ----------------------------------------------------------------------------- 
// Author:  James Lambert
// Created: February 2nd, 2018
// 
// Class to measure the configuration during the simulation and output the
// results to a file for processing.
// -----------------------------------------------------------------------------

#ifndef PROBE_H
#define PROBE_H

#include <vector>
#include <string>
#include "CONFIG.h"
#include "OBSERVABLE.h"

namespace SSE
{
  class PROBE
  {
    private:
      
      int     _ns; 
      int     _nb;
      double  _df;
      double  _bt;
      double  _ep;
      bool    _bc;
      double  _offset;
      double  _kval;
     

      // These switches control which quantities are measured by default
      bool _menergy    =   false;
      bool _msqarstrg  =   false;
      bool _msqarspin  =   false;
      bool _mcorrfunc  =   false;
      bool _msqarcorr  =   false;
      bool _mspinavrg  =   false;
      bool _mstrgcorr  =   false;

      OBSERVABLE<double>               _energy;
      OBSERVABLE<double>               _nave;
      OBSERVABLE<double>               _nsquared;
      OBSERVABLE<double>               _sqarspin;
      OBSERVABLE<double>               _smsusc;
      OBSERVABLE<double>               _sosusc;
      OBSERVABLE<double>               _sm;
      OBSERVABLE<double>               _smsq;
      OBSERVABLE<double>               _so;
      OBSERVABLE<double>               _sosq;
      std::vector<OBSERVABLE<double> > _corrfunc;

      int sfac(int i){
        if(i%2 == 0) return 1;
        else return -1;
      }

    public:

      PROBE(CONFIG& c, 
            const double& k
           ) : _ns(c.ns()), 
               _nb(c.nb()), 
               _df(c.df()), 
               _bt(c.bt()), 
               _ep(c.ep()),  
               _bc(c.bc()), 
               _offset(c.eo()), 
               _kval(k) {}
      
      PROBE(CONFIG& c,        
            const double& k,  // value of momentum for spin operator  
            const bool& A,    // energy switch
            const bool& B,    // string correlation squared switch
            const bool& C,    // average of square of sum of spins
            const bool& D,    // Sz correlations switch
            const bool& E,    // squared Sz spin averages 
            const bool& F,    // spin averages
            const bool& G     // string correlations switch
           ) : _ns(c.ns()), 
               _nb(c.nb()), 
               _df(c.df()),
               _bt(c.bt()), 
               _ep(c.ep()), 
               _bc(c.bc()),
               _offset(c.eo()),
               _kval(k), 
               _menergy(A),
               _msqarstrg(B),
               _msqarspin(C), 
               _mcorrfunc(D), 
               _msqarcorr(E),
               _mspinavrg(F), 
               _mstrgcorr(G){}
     
      void meas_energy(CONFIG&);
      void meas_nave(CONFIG&);
      void meas_nsquared(CONFIG&);
      void meas_qvar(CONFIG&);
      void meas_corrfunc(CONFIG&);

      void output_to_file(const std::string&);
      
      void operator()(CONFIG& c)
      {
        if(_menergy)   meas_energy(c);
        meas_nave(c);
        meas_nsquared(c);
        if(_mcorrfunc) meas_corrfunc(c);
      }
  };
}

#endif // PROBE_H
