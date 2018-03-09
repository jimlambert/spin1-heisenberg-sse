// ----------------------------------------------------------------------------- 
// Author:  James Lambert
// Created: January 28th, 2018
//
// Class which contains all key data structures for the stochastic series
// expansion for the spin-1 AFM Heisenberg model with uniaxial anisotropy, as 
// well as the update routines.
// -----------------------------------------------------------------------------

#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <array>
#include <random>

namespace SSE
{
  class CONFIG
  {
    private:
      
      int               _ns;        // number of sites in lattice
      int               _nb;        // number of bonds in lattice
      int               _xo;        // current expansion order of lattice
      int               _no;        // current number of operators in lattice
      double            _bt;        // inverse temperature in units of k_b
      double            _df;        // single ion anisotropy term
      double            _ep;        // small parameter to add to weights
      bool              _bc;        // boundary conditions, PBC=true
      int               _pi;        // current place in the propagation list

      int*              _spins;     // contains spins on each bond site
      int*              _sites[2];  // contains structure of lattice
      std::vector<int>  _oplst;     // bonds on which operators act
      std::vector<int>  _vtlst;     // list of vertex types
      double            _vwgts[17]; // vertex weights needed for updates

      // lookup tables for exit vertices and exit probability bounds for the
      // single spin flip loop update
      int               _outvrts[136][4];
      double            _extprbs[136][4];

      // look up tables for exit vertices and exit probability boudns for the
      // double spin flip loop update
      int               _doutvrts[68][4];
      double            _dextprbs[68][4]; 
      double            _dbfprbs[2];

      // seed random number generator and initialize distributions
      std::random_device _rd;
      std::mt19937 _mteng{_rd()};
      std::uniform_real_distribution<> _rdist{0.0, 0.99999};
      std::uniform_int_distribution<>  _rspin{-1, 1};
      std::uniform_int_distribution<>  _rleg{0, 3};
      std::uniform_int_distribution<>  _rud{0, 1};
      
      // needs to be initialized after number of bonds is known
      std::uniform_int_distribution<>* _rbond;

      // function determines the probabily index given an input leg, an input
      // vertex, and a change type. Function is ideally inlined
      int _prbindex(const int& e, const int& v, const int& ds)
      {
        return e + ((v - 1) * 8) + 4 * ds;
      }

      // same as above but now for the double spin flips
      int _dprbindex(const int& e, const int& v)
      {
        return e + ((v - 1) * 4);
      }

      // given four spins, returns the corresponding vertex id. Also inlined
      // ideally due to number of function calls. If an illegal vertex is
      // accessed an id of 0 is return which is likely cause a seg fault if not
      // handled correctly.
      int _vrtid(const int (&spins)[4]) const
      {
        int s1 = spins[0], s2 = spins[1], s3 = spins[2], s4 = spins[3];

        if     ((s1 == 0) && (s2 == 0) && (s3 == 0) && (s4 == 0)) return 1;
        else if((s1 == 1) && (s2 == 0) && (s3 == 1) && (s4 == 0)) return 2;
        else if((s1 ==-1) && (s2 == 0) && (s3 ==-1) && (s4 == 0)) return 3;
        else if((s1 == 0) && (s2 == 1) && (s3 == 0) && (s4 == 1)) return 4;
        else if((s1 == 0) && (s2 ==-1) && (s3 == 0) && (s4 ==-1)) return 5;
        else if((s1 == 1) && (s2 == 1) && (s3 == 1) && (s4 == 1)) return 6;
        else if((s1 ==-1) && (s2 ==-1) && (s3 ==-1) && (s4 ==-1)) return 7;
        else if((s1 == 1) && (s2 ==-1) && (s3 == 1) && (s4 ==-1)) return 8;
        else if((s1 ==-1) && (s2 == 1) && (s3 ==-1) && (s4 == 1)) return 9;
        else if((s1 == 0) && (s2 == 0) && (s3 == 1) && (s4 ==-1)) return 10;
        else if((s1 == 0) && (s2 == 0) && (s3 ==-1) && (s4 == 1)) return 11;
        else if((s1 ==-1) && (s2 == 1) && (s3 == 0) && (s4 == 0)) return 12;
        else if((s1 == 1) && (s2 ==-1) && (s3 == 0) && (s4 == 0)) return 13;
        else if((s1 == 1) && (s2 == 0) && (s3 == 0) && (s4 == 1)) return 14;
        else if((s1 ==-1) && (s2 == 0) && (s3 == 0) && (s4 ==-1)) return 15;
        else if((s1 == 0) && (s2 == 1) && (s3 == 1) && (s4 == 0)) return 16;
        else if((s1 == 0) && (s2 ==-1) && (s3 ==-1) && (s4 == 0)) return 17;
        else                                                      return 0;
      }

      // This function will only return diagonal vertices. It's used in the
      // diagonal update which is called many thousands of times and should be
      // faster than the above funcion.
      int _diagvrtid(const int& s1, const int& s2)
      {
        if     ((s1 == 0) && (s2 == 0)) return 1;
        else if((s1 == 1) && (s2 == 0)) return 2;
        else if((s1 ==-1) && (s2 == 0)) return 3;
        else if((s1 == 0) && (s2 == 1)) return 4;
        else if((s1 == 0) && (s2 ==-1)) return 5;
        else if((s1 == 1) && (s2 == 1)) return 6;
        else if((s1 ==-1) && (s2 ==-1)) return 7;
        else if((s1 == 1) && (s2 ==-1)) return 8;
        else                            return 9;
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
      
      ~CONFIG();

      void expoupdt();     // expansion order update
      void diagupdt();     // diagonal operator update
      void loopupdt();     // directed loop update

      // propagates internal spins state according to the operator list. Useful
      // for making measurements of the correlations functions and averaging
      // over the configuration space.
      void propagate(); 

      // access operator for spin chain
      int operator[](const int i) const{return _spins[i];}
      
      // access members for internal configuration variables
      int ns()      const {return _ns;}
      int nb()      const {return _nb;}
      int xo()      const {return _xo;}
      int no()      const {return _no;}
      double bt()   const {return _bt;}
      double df()   const {return _df;}
      double ep()   const {return _ep;}
      bool bc()     const {return _bc;}
      int op(int p) const {return _vtlst[p];}
      double eo()   const {return 1 + _df + _ep;}

      // display functions - useful for debugging
      void disp_wgts();
      void disp_outvrts(const bool&);
      void disp_extprbs(); 
      void disp_spins();
      void disp_config();
      void disp_opers();
      void disp_vert(const int&);
      
      // the configuration arrays are hard coded and carried around by the
      // object as public members so that the measurements functions may use
      // them if necessary 
      const int verts[17][4] = 
      {
        { 0,   0,  0,  0},  // 1
        { 1,   0,  1,  0},  // 2
        {-1,   0, -1,  0},  // 3
        { 0,   1,  0,  1},  // 4 
        { 0,  -1,  0, -1},  // 5 
        { 1,   1,  1,  1},  // 6
        {-1,  -1, -1, -1},  // 7 
        { 1,  -1,  1, -1},  // 8
        {-1,   1, -1,  1},  // 9
        { 0,   0,  1, -1},  // 10
        { 0,   0, -1,  1},  // 11 
        {-1,   1,  0,  0},  // 12
        { 1,  -1,  0,  0},  // 13
        { 1,   0,  0,  1},  // 14
        {-1,   0,  0, -1},  // 15
        { 0,   1,  1,  0},  // 16 
        { 0,  -1, -1,  0}   // 17
      };

      const bool types[17] = 
      {
        true,               // 1
        true,               // 2
        true,               // 3
        true,               // 4
        true,               // 5
        true,               // 6
        true,               // 7
        true,               // 8
        true,               // 9
        false,              // 10 
        false,              // 11
        false,              // 12
        false,              // 13
        false,              // 14
        false,              // 15
        false,              // 16
        false               // 17
      };

      // Output flips for the completely deterministic two spin flip update.
      const int dbouts[17][2] =
      {
        {0, 0},             // 1
        {3, 0},             // 2
        {2, 0},             // 3
        {0, 5},             // 4
        {0, 4},             // 5
        {9, 8},             // 6
        {8, 9},             // 7
        {7, 6},             // 8
        {6, 7},             // 9
        {11, 11},           // 10
        {10, 10},           // 11
        {13, 13},           // 12
        {12, 12},           // 13
        {15, 15},           // 14
        {14, 14},           // 15
        {17, 17},           // 16
        {16, 16}            // 17  
      };

      const int dbexts[3][4] =
      {
        {2, 3, 0, 1},   
        {1, 0, 3, 2},
        {3, 2, 1, 0}
      }; 

  };
}

#endif  // CONFIG_H
