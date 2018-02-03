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
#include "CONFIG.h"
#include "observable.h"

namespace SSE
{
  class PROBE
  {
    private:

      int    _ns; 
      int    _nb;
      double _df;
      double _bt;
      double _ep;

      CONFIG* conf;

      SSE::observable<double>               _energy;
      std::vector<SSE::observable<double> > _corrfunc;

    public:

      PROBE();

      void meas_energy();
  };
}

#endif // PROBE_H
