#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <array>
#include "CONFIG.h"

// VRTTABLE   -   allowed vertices for the spin-1 AFM Heisenberg model with
//                uniaxial anisotropy
// TYPTABLE   -   Contains array that indicates whether the vertex at the 
//                corresponding index is diagonal or off-diagonal
// OUTTABLE   -   contains output vertices each corresponding to exiting on a
//                given leg of the bare vertex
// PRBTABLE   -   contains the bounds for the transition probabilities given an
//                input leg, input vertex, and change type.

//typedef std::array<<std::array<int, 4>,    17>   VRTTABLE;   
//typedef std::array<bool,                   17>   TYPTABLE;
//typedef std::array<<std::array<int, 4>,    136>  OUTTABLE;   
//typedef std::array<<std::array<double, 4>, 136>  PRBTABLE;   

//std::array<<std::array<int, 4>, 17> verts = 
static const int verts[17][4] = 
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

static const bool types[17] = 
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


// define display functions for the above structures (these could be folded into
// the class... oh well.
void disp_vert(const int& j)
{
  int i = j-1;
  if(j==0) return;
  std::cout << "vertex id:" << '\t' << j << std::endl;
  std::cout << std::setw(5) << std::setfill(' ') 
            << std::left << verts[i][0] << '\t' << verts[i][1] 
            << std::endl 
            << std::setw(10) << std::setfill('-') << '-' << std::endl 
            << std::setw(5) << std::setfill(' ') 
            << std::left << verts[i][2] << '\t' << verts[i][3] 
            << std::endl << std::endl;
}

namespace SSE
{
  CONFIG::CONFIG(const int& N, const double& T, const double& D, 
                 const double& E, const bool& BC) : 
                 _ns(N), _xo(20), _no(0), _bt(1.0/T), _df(D), _ep(E), _bc(BC)
  {
    // seed random number generator
    srand(time(NULL));

    // initialize spin chain with random spins
    _spins = new int[_ns];
    for(int i=0; i<_ns; i++)
    {
      double r = (double)rand() / (double)RAND_MAX;
      if     (r < 0.33) _spins[i] = -1;
      else if(r < 0.67) _spins[i] = 0;
      else              _spins[i] = 1;
    }
    
    // initialize the sites array which contains the bond structure
    if(_bc)   // PBC
    {
      _nb = _ns;
      _sites[0] = new int[_nb];
      _sites[1] = new int[_nb];

      for(int i=0; i<_nb; i++)
      {
        _sites[0][i] = i;
        if(i<_nb - 1) _sites[1][i] = i+1;
        else          _sites[1][i] = 0;
      } 
    }
    else      // OBC
    {
      _nb = _ns - 1;
      _sites[0] = new int[_nb];
      _sites[1] = new int[_nb];

      for(int i=0; i<_nb; i++)
      {
        _sites[0][i] = i;
        _sites[1][i] = i + 1;
      }
    }
    
    // start with an initial expansion order of 20. In both _oplist and _vtlist
    // the number 0 denotes the identity vertex.
    _oplst.resize(_xo, 0);
    _vtlst.resize(_xo, 0);

    // we assume J=1. Using this, that constant offset guarenteeing positive
    // definiteness is given by 1 + D + eps. We assume a hamiltonian with bond
    // operator H_b = J S_i^zS_j^z + D/2((S_i^z)^2 + (S_j^z)^2) + O.D. where
    // O.D. represents the off diagonal term. 
    double C = 1 + _df + _ep;
      
    for(unsigned int i=0; i<17; i++)
    {
      if(types[i])                      // diagonal
      {
        int s1=verts[i][0], s2 = verts[i][1];
        _vwgts[i] = C - (s1*s2) - 0.5 * _df * (s1*s1 + s2*s2);
      }
      else _vwgts[i] = 1.0;            // off-diagonal    
    } 
    
    // populate look-up table for transition probabilities
    _calcprobs(); 
  }
  
  void CONFIG::_calcprobs()
  {
    // first determine all output vertices
    for(unsigned int ref=0; ref<17; ref++)
    {
      // loop through all in and out legs starting with up flips on the entrance
      // followed by down flips on the entrance
      for(unsigned int i=0; i<4; i++){
        for(unsigned int o=0; o<4; o++)
        {
          int vspu[4] = {verts[ref][0], verts[ref][1], 
                         verts[ref][2], verts[ref][3]};
          int vspd[4] = {verts[ref][0], verts[ref][1], 
                         verts[ref][2], verts[ref][3]};
          vspu[i] = vspu[i] + 1;
          if((i<2 && o<2) || (i>1 && o>1)) vspu[o] = vspu[o] - 1;
          else                             vspu[o] = vspu[o] + 1;
          _outvrts[_prbindex(i, ref+1, 1)][o] = _vrtid(vspu);
          vspd[i] = vspd[i] - 1;
          if((i<2 && o<2) || (i>1 && o>1)) vspd[o] = vspd[o] + 1;
          else                             vspd[o] = vspd[o] - 1;
          _outvrts[_prbindex(i, ref+1, -1)][o] = _vrtid(vspd);
        }
      } 
    }
    // now we determine transition probabilities and bounds
    for(unsigned int row=0; row<136; row++)
    {
      int x;                        // exit leg index
      double den = 0.0;             // denominator for heatbath calculation  
      double sum = 0.0;             // used to calculate probability bounds
      std::array<double, 4> temp;   // temporary array to store porbabilies
      for(x=0; x<4; x++) den += _vwgts[_outvrts[row][x]-1];      
      for(x=0; x<4; x++)
      {
        if(_outvrts[row][x] == 0) temp[x] = 0.0;
        else temp[x] = _vwgts[_outvrts[row][x]-1] / den;
      }

      int lstval = 0;
      // convert these probabilities to bounds
      for(x=0; x<4; x++)
      {
        sum += temp[x];
        if(temp[x] == 0.0) _extprbs[row][x] = 0.0;
        else
        {
          lstval = x;
          _extprbs[row][x] = sum;
        }         
        // guarentee that last non-zero bound is 1. 
      }  
      _extprbs[row][lstval] = 1.0;
    }
  }
 
  void CONFIG::expoupdt()
  {
    // called during the equilibration phase. This update increases the total
    // possible expansion order of the system.
    int newxo = 1.33 * _no;
    if(newxo > _xo)
    {
      _oplst.resize(newxo, 0);
      _vtlst.resize(newxo, 0);
      _xo = newxo;
    }
  }

  void CONFIG::diagupdt()
  {
    // loop through the configuration and decide to insert of remove an operator
    for(int p=0; p<_xo; p++)
    {
      if(_oplst[p] == 0)              // no operator present
      {
        // select a random bond to attempt insertion
        int rand_bond = rand() % _nb;
        int spin1 = _spins[_sites[0][rand_bond]];
        int spin2 = _spins[_sites[1][rand_bond]];
        int id = _diagvrtid(spin1, spin2);

        // compare random number to weight
        double r = (double)rand() / (double)RAND_MAX;
        if((r*(_xo-_no)<_nb*_bt*_vwgts[id-1]))
        {
          _oplst[p] = (1 + rand_bond) * 2;
          _vtlst[p] = id;
          _no += 1;
        }

      }
      else if(types[_vtlst[p]-1])     // diagonal operator present
      {
        // compare random number to removal probability
        double r = (double)rand() / (double)RAND_MAX; 
        int id = _vtlst[p];
        if((r*_nb*_bt*_vwgts[id-1])<(_xo-_no+1))
        {
          _oplst[p] = 0;
          _vtlst[p] = 0;
          _no -= 1;
        }
      }
      else                            //  off-diagonal operator present
      {
        // propagate the internal spin state
        int bond = int(_oplst[p] / 2) - 1;
        int spin1 = verts[_vtlst[p]][2];
        int spin2 = verts[_vtlst[p]][2];
        _spins[_sites[0][bond]] = spin1;
        _spins[_sites[1][bond]] = spin2; 
      }
    }
  }

  void CONFIG::disp_wgts()
  {
    for(int i=0; i<17; i++)
    {
      std::cout << "vertex: " << i + 1 << " -> " << _vwgts[i] << std::endl;
    }
  }
  
  void CONFIG::disp_outvrts(const bool& image)
  {
    if(image)
    {
      for(int i=0; i<136; i++)
      {
        std::cout << std::setw(20) << std::setfill('-') << "-" << std::endl;    
        if(i % 8 > 4) std::cout << "FLIP UP" << std::endl;
        else std::cout << "FLIP DOWN" << std::endl;
        std::cout << "entrance leg: " << i % 4 << std::endl;
        disp_vert(i/8 + 1);
        std::cout << std::setw(10) << std::setfill('=') << "=" << std::endl;
        for(int j=0; j<4; j++)
        {
          if(_outvrts[i][j] == 0) continue;
          std::cout << "exit leg: " << j << std::endl;
          disp_vert(_outvrts[i][j]);
          std::cout << std::endl;
        }
      } 
    }
    else
    {
      for(int i=0; i<136; i++){
        for(int j=0; j<4; j++){
          std::cout << _outvrts[i][j] <<'\t';
        }
      }
    }  
  }

  void CONFIG::disp_extprbs()
  {
    for(int i=0; i<136; i++){
      for(int j=0; j<4; j++){
        std::cout << _extprbs[i][j] << '\t';
      }
      std::cout << std::endl;
    }
  }

  void CONFIG::disp_spins()
  {
    for(int i=0; i<_ns; i++)
    {
      if     (_spins[i] ==-1) std::cout << "-" << '\t';
      else if(_spins[i] == 0) std::cout << "0" << '\t';
      else                    std::cout << "+" << '\t';
    }
    std::cout << std::endl;
  }

  void CONFIG::disp_config()
  {
    disp_spins();
    for(int p=0; p<_xo; p++)
    {
      if(_oplst[p] == 0)
      {
        std::cout << std::endl;
        disp_spins();
        continue;
      }
      int bond = (_oplst[p] / 2) - 1;
      for(int i=0; i<bond; i++)
      {
        std::cout << '\t';
      }
      if(types[_vtlst[p]-1]) std::cout << "000000000" << std::endl;
      else std::cout << "XXXXXXXXX" << std::endl;
      disp_spins();
    }
  }

  void CONFIG::disp_opers()
  {
    for(int p=0; p<_xo; p++)
    {
      std::cout << "[" << p << "]" 
                << '\t' << _oplst[p] 
                << '\t' << _vtlst[p] 
                << std::endl;
    }
  }

}
