#include <iostream>
#include <iomanip>
#include <array>
#include <vector>
#include <random>
#include "CONFIG.h"

namespace SSE
{
  CONFIG::CONFIG(const int& N, const double& T, const double& D, 
                 const double& E, const bool& BC) : 
                 _ns(N), _xo(20), _no(0), _bt(1.0/T), _df(D), _ep(E), _bc(BC),
                 _pi(0)
  {
    // initialize spin chain with random spins
    _spins = new int[_ns];
    for(int i=0; i<_ns; i++)
      _spins[i] = _rspin(_mteng); 
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

    _rbond = new std::uniform_int_distribution<>{0, _nb-1}; 
  
    // start with an initial expansion order of 20. In both _oplist and _vtlist
    // the number 0 denotes the identity vertex.
    _oplst.resize(_xo, 0);
    _vtlst.resize(_xo, 0);

    // we assume J=1. Using this, that constant offset guarenteeing positive
    // definiteness is given by 1 + D + eps. We assume a hamiltonian with bond
    // operator H_b = J S_i^zS_j^z + D/2((S_i^z)^2 + (S_j^z)^2) + O.D. where
    // O.D. represents the off diagonal term. 
    double C = 1 + _df + _ep;
      
    for(unsigned int i=0; i<17; i++){
      if(types[i])                      // diagonal
      {
        int s1=verts[i][0], s2 = verts[i][1];
        _vwgts[i] = C - (s1*s2) - 0.5 * _df * (s1*s1 + s2*s2);
      }
      else _vwgts[i] = 1.0;            // off-diagonal    
    } 
    
    // populate look-up table for transition probabilities
    _calcprobs();
    // for the double flip update we need the transition probability between the
    // ferromagnetic and antiferromagnetic vertices. This is achieved with a
    // heat bath algorithm.
    double tot = _vwgts[6] + _vwgts[7];
    _dbfprbs[0] = _vwgts[6] / tot;
    _dbfprbs[1] = _vwgts[7] / tot;
  }
  
  CONFIG::~CONFIG()
  {
    delete [] _spins;
    delete [] _sites[0];
    delete [] _sites[1];
    delete    _rbond;
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
          // single flip update ------------------------------------------------
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
          _outvrts[_prbindex(i, ref+1, 0)][o] = _vrtid(vspd);

          // double flip update ------------------------------------------------
          int vsp[4] = {verts[ref][0], verts[ref][1], 
                        verts[ref][2], verts[ref][3]};
          if     (vsp[o] == 1)  vsp[o] = vsp[o] - 2;
          else if(vsp[o] == -1) vsp[o] = vsp[o] + 2;
          if     (vsp[i] == 1)  vsp[i] = vsp[i] - 2;
          else if(vsp[i] == -1) vsp[i] = vsp[i] + 2;
          _doutvrts[_dprbindex(i, ref+1)][o] = _vrtid(vsp);
        }
      }
    }
    // single flip update ------------------------------------------------------
    for(unsigned int row=0; row<136; row++)
    {
      int x;                        // exit leg index
      double den = 0.0;             // denominator for heatbath calculation  
      double sum = 0.0;             // used to calculate probability bounds
      std::vector<double> temp;     // temporary array to store probabilities
      temp.resize(4, 0.0);
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
      }  
      _extprbs[row][lstval] = 1.0;
    }
    
    // double flip update ------------------------------------------------------
    for(unsigned int row=0; row<68; row++)
    {
      int x;                        // exit leg index
      double den = 0.0;             // denominator for heatbath calculation  
      double sum = 0.0;             // used to calculate probability bounds
      std::vector<double> temp;     // temporary array to store probabilities
      temp.resize(4, 0.0);
      for(x=0; x<4; x++) den += _vwgts[_doutvrts[row][x]-1];
      for(x=0; x<4; x++)
      {
        if(_doutvrts[row][x] == 0) temp[x] = 0.0;
        else temp[x] = _vwgts[_doutvrts[row][x]-1] / den;
      }
      int lstval = 0;
      // convert probabilities to bounds
      for(x=0; x<4; x++)
      {
        sum += temp[x];
        if(temp[x] == 0.0) _dextprbs[row][x] = 0.0;
        else
        {
          lstval = x;
          _dextprbs[row][x] = sum;
        }
      }
      _dextprbs[row][lstval] = 1.0; 
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
      if(_vtlst[p] == 0)              // no operator present
      {
        // select a random bond to attempt insertion
        int rand_bond = (*_rbond)(_mteng);
        int spin1 = _spins[_sites[0][rand_bond]];
        int spin2 = _spins[_sites[1][rand_bond]];
        int id = _diagvrtid(spin1, spin2);

        // compare random number to weight
        double r = _rdist(_mteng);
        if((r*(_xo-_no)<_nb*_bt*_vwgts[id-1]))
        {
          _oplst[p] = rand_bond;
          _vtlst[p] = id;
          _no += 1;
        }

      }
      else if(types[_vtlst[p]-1])     // diagonal operator present
      {
        // compare random number to removal probability
        double r = _rdist(_mteng); 
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
        int bond = _oplst[p];
        int spin1 = verts[_vtlst[p]-1][2];
        int spin2 = verts[_vtlst[p]-1][3];
        _spins[_sites[0][bond]] = spin1;
        _spins[_sites[1][bond]] = spin2; 
      }
    }
  }

  void CONFIG::loopupdt()
  {
    // initialize list containing information about vertex links
    std::vector<int> linklst; 
    std::vector<int> vfrst;
    std::vector<int> vlast;
    linklst.resize(4*_xo+1, 0);
    vfrst.resize(_ns, -1);
    vlast.resize(_ns, -1);

    int v0, v1, v2, s1, s2;

    for(int p=0; p<_xo; p++)
    {
      v0 = 4 * p + 1;
      if(_vtlst[p] == 0){
        for(int i=0; i<4; i++) linklst[v0 + i] = 0;
      }
      else
      { 
        s1 = _sites[0][_oplst[p]];
        s2 = _sites[1][_oplst[p]];
        v1 = vlast[s1];
        v2 = vlast[s2];
        if(v1 != -1)
        {
          linklst[v1] = v0; 
          linklst[v0] = v1;
        }
        else vfrst[s1] = v0;
        if(v2 != -1)
        {
          linklst[v2]   = v0 + 1;
          linklst[v0+1] = v2;
        }
        else vfrst[s2] = v0 + 1;
        vlast[s1] = v0 + 2;
        vlast[s2] = v0 + 3;
      }
    }

    // close links across imaginary time boundary
    for(int i=0; i<_ns; i++){
      if(vfrst[i] != -1){
        linklst[vlast[i]] = vfrst[i];
        linklst[vfrst[i]] = vlast[i];
      }
    }
    
    /*  Print vertex list for debugging purposes
    for(int i=1; i<(4*_xo+1); i++)
    {
      std::cout << "[" << i << "]" << '\t' << linklst[i]
                << '\t';
      if(i%4 == 0) std::cout << std::endl;
    }
    */
    
    // create temporary list to store changes to vertex list 
    std::vector<int> newvrts;
    for(unsigned int i=0; i<_vtlst.size(); i++)
      newvrts.push_back(_vtlst[i]);

    for(int p=0; p<_xo; p++)
    {
      if(_vtlst[p] == 0) continue;
      else
      {
        // Choose starting leg
        int e = _rleg(_mteng);
        int v0 = 4 * p + 1 + e;
        int vc = v0;
        int ud, ut = _rud(_mteng);
        // Determine the initial direction of the spin flip
        if(verts[newvrts[p]-1][e] == 0)      ud = _rud(_mteng);
        else if(verts[newvrts[p]-1][e] == 1) ud = 0;
        else                                 ud = 1;

        //if(abs(verts[newvrts[p]-1][e])==1 && (0.0 > _rdist(_mteng))) continue;
        do
        { 
          // generate a random number here
          double r = _rdist(_mteng);
          int x=e;
          int ind = _prbindex(e, newvrts[(vc-1)/4], ud);
          // choose exit leg
          for(int i=0; i<4; i++){
            if(r<_extprbs[ind][i]){
              x=i; 
              break;
            }
          }
          // determine if spin flip is up or down.
          if((e+x==5)||(e+x==1)) ud = ud^1;
          if(x==e)               ud = ud^1;
          // change vertex 
          newvrts[(vc-1)/4] = _outvrts[ind][x];
          // move to new leg
          vc = linklst[vc - e + x];
          e = (vc-1) % 4;
        }while(vc != v0);
      }
    }
    // commit changes
    for(unsigned int i=0; i<_vtlst.size(); i++)
      _vtlst[i] = newvrts[i];

    // flip spins without operators
    for(int i=0; i<_ns; i++)
    {
      if(vfrst[i] != -1)
      {
         int index = (vfrst[i]-1) / 4;
         int leg   = (vfrst[i]-1) % 4;
         _spins[i] = verts[_vtlst[index]-1][leg];
      }
      else _spins[i] = _rspin(_mteng);    // generate new spin
    }
  }
 
  void CONFIG::propagate()
  {
    if(_vtlst[_pi] == 0)
    {
      _pi = (_pi + 1) % _xo;
      return;
    }
    if(!types[_vtlst[_pi]-1]) 
    {
      _spins[_sites[0][_oplst[_pi]]] = verts[_vtlst[_pi]-1][2];
      _spins[_sites[1][_oplst[_pi]]] = verts[_vtlst[_pi]-1][3];
    }
    _pi = (_pi + 1) % _xo;
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
        std::cout << "[" << i << "][";
        for(int j=0; j<4; j++){
          std::cout << _outvrts[i][j] << " ";
        }
        std::cout << "]" << std::endl;
      }
    }  
  }

  void CONFIG::disp_extprbs()
  {
    for(int i=0; i<136; i++){
      int vid = (i / 8 + 1);
      if(i%8==0) std::cout << "vertex: " << vid << std::endl;
      for(int j=0; j<4; j++){
        std::cout << _extprbs[i][j] << '\t';
      }
      std::cout << std::endl;
      if(i%8==7) std::cout << "------------------------" << std::endl;
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
      if(_vtlst[p] == 0)
      {
        std::cout << std::endl;
        disp_spins();
        continue;
      }
      int bond = _oplst[p];
      for(int i=0; i<bond; i++)
      {
        std::cout << '\t';
      }
      if(types[_vtlst[p]-1]) std::cout << "000000000";
      else
      { 
        std::cout << "XXXXXXXXX";
        _spins[_sites[0][_oplst[p]]] = verts[_vtlst[p]-1][2];
        _spins[_sites[1][_oplst[p]]] = verts[_vtlst[p]-1][3];
      } 
      for(int i=0; i<(_nb - bond + 1); i++)
      {
        std::cout << '\t';
      }
      std::cout << "[" << p << "] " << _vtlst[p] << std::endl;
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
  
  void CONFIG::disp_vert(const int& j)
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
}
