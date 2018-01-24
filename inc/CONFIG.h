// ----------------------------------------------------------------------------- 
// Author:  James Lambert
// Created: Januaray 28th, 2018
//
// Class which contains all key data structures for the stochastic series
// expansion for the spin-1 AFM Heisenberg model with uniaxial anisotropy, as 
// well as the update routines.
// -----------------------------------------------------------------------------

#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <array>

// These functions are used to display global source file variables
void disp_vert(const int&);
void disp_outvrts();
void disp_extprbs();


namespace SSE
{
  class CONFIG
  {
    private:
      
      int     _ns;                    // number of sites in lattice
      int     _nb;                    // number of bonds in lattice
      int     _xo;                    // current expansion order of lattice
      int     _no;                    // current number of operators in lattice
      double  _bt;                    // inverse temperature in units of k_b
      double  _df;                    // single ion anisotropy term
      double  _ep;                    // small parameter to add to weights
      bool    _bc;                    // boundary conditions, PBC=true

      int*                _spins;     // contains spins on each bond site
      int*                _sites[2];  // contains structure of lattice
      std::vector<int>    _oplst;     // bonds on which operators act
      std::vector<int>    _vtlst;     // list of vertex types
      double              _vwgts[17]; // vertex weights needed for updates


      // function determines the probabily index given an input leg, an input
      // vertex, and a change type. Function is ideally inlined
      int _prbindex(const int& e, const int& v, const int& ds)
      {
        return e + (v - 1) * 8 + int(ds + 1.5) * 4;
      }

      // given four spins, returns the corresponding vertex id. Also inlined
      // ideally due to number of function calls. If an illegal vertex is
      // accessed an id of 0 is return which is likely cause a seg fault if not
      // handled correctly.
      int _vrtid(const int (&spins)[4]) const
      {
        int s1 = spins[0], s2 = spins[1], s3 = spins[2];

        if     ((s1 == 0) && (s2 == 0) && (s3 == 0)) return 1;
        else if((s1 == 1) && (s2 == 0) && (s3 == 1)) return 2;
        else if((s1 ==-1) && (s2 == 0) && (s3 ==-1)) return 3;
        else if((s1 == 0) && (s2 == 1) && (s3 == 0)) return 4;
        else if((s1 == 0) && (s2 ==-1) && (s3 == 0)) return 5;
        else if((s1 == 1) && (s2 == 1))              return 6;
        else if((s1 ==-1) && (s2 ==-1))              return 7;
        else if((s1 == 1) && (s2 ==-1) && (s3 == 1)) return 8;
        else if((s1 ==-1) && (s2 == 1) && (s3 ==-1)) return 9;
        else if((s1 == 0) && (s2 == 0) && (s3 == 1)) return 10;
        else if((s1 == 0) && (s2 == 0) && (s3 ==-1)) return 11;
        else if((s1 ==-1) && (s2 == 1) && (s3 == 0)) return 12;
        else if((s1 == 1) && (s2 ==-1) && (s3 == 0)) return 13;
        else if((s1 == 1) && (s2 == 0) && (s3 == 0)) return 14;
        else if((s1 ==-1) && (s2 == 0) && (s3 == 0)) return 15;
        else if((s1 == 0) && (s2 == 1) && (s3 == 1)) return 16;
        else if((s1 == 0) && (s2 ==-1) && (s3 ==-1)) return 17;
        else                                         return 0;
      }     

      // calculate the transition probabilities for the update.
      void _calcprobs();
      
    public:
 
      CONFIG(const int&,      // system size
             const double&,   // temperature
             const double&,   // D-field strength
             const double&,   // value of epsilon
             const bool&      // true for periodic boudary conditions
            );
      
      void expo_update();     // expansion order update
      void diag_update();     // diagonal operator update
      void loop_update();     // directed loop update

      // propagates internal spins state according to the operator list. Useful
      // for making measurements of the correlations functions and averaging
      // over the configuration space.
      void propagate(); 

      // access operator for spins
      int operator[](const int& i) const{return _spins[i];}
      int ns() const {return _ns;}
      int xo() const {return _xo;}
      int no() const {return _no;} 
  };
}

#endif  // CONFIG_H
