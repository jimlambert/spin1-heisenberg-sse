#include <vector>
#include <complex>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "OBSERVABLE.h"
#include "CONFIG.h"
#include "PROBE.h"

#define CORR_INT 1

namespace SSE
{
  void PROBE::meas_energy(CONFIG& conf)
  {
    double val = (-1.0/_bt)*(conf.no()) + _ns*_offset;
    _energy.push_back(val);
  }
  
  void PROBE::meas_nave(CONFIG& c){
    _nave.push_back(c.no());
  }

  void PROBE::meas_nsquared(CONFIG& c){
    _nsquared.push_back(c.no()*c.no());
  }
  
  void PROBE::meas_qvar(CONFIG& c){
    double smsc = 0.0;
    double smsq = 0.0;
    double sosc = 0.0;
    double sosq = 0.0;
    int meascount = 0;
    for(int p=0; p<c.xo(); p++){
      if(p % CORR_INT == 0){
        double sctot = 0.0;
        double sotot = 0.0;
        int sptot = 0.0;
        for(int i=0; i<c.ns(); i++){
          sctot += 0.5*sfac(i)*c[i];
          sotot += c[0]*c[i]*sfac(sptot);
          sptot += c[i];
        }
        smsc += sctot;
        smsq += sctot*sctot;
        sosc += sotot;
        sosq += sotot*sotot; 
        meascount += 1;
      }
      c.propagate();
    }
    _sm.push_back(smsc/meascount);
    _so.push_back(sosc/meascount);
    double smval = c.bt()/(meascount*(meascount+1))*smsc*smsc 
                 + c.bt()/(meascount*(meascount+1))*smsq;
    _smsusc.push_back(smval-c.bt()*(smsc*smsc/(meascount*meascount)));
    _smfluc.push_back(smsq/(meascount) - (smsc*smsc)/((meascount)*(meascount)));   
    double soval = c.bt()/(meascount*(meascount+1))*sosc*sosc
                 + c.bt()/(meascount*(meascount+1))*sosq;
    _sosusc.push_back(soval-c.bt()*(sosc*sosc/(meascount*meascount)));
    _sofluc.push_back(sosq/meascount - (sosc*sosc)/(meascount*meascount));
  }

  void PROBE::meas_corrfunc(CONFIG& conf)
  {
    int max;
    int initial;
    if(_bc){  // PBC
      initial = 0;
      max = _ns / 2 + 1;
    } 
    else{     // OBC
      initial = _ns /2;
      max = _ns;
    }    
    // if the corr func vector is of size zero it needs to be initialized with
    // the first observables
    if(_corrfunc.empty())
    {
      for(int i=initial; i<max; i++)
      {
        OBSERVABLE<double> tempobs;
        double total = 0.0;
        for(int p=0; p<(conf.xo()); p++)
        {
          total += conf[0] * conf[i];
          conf.propagate();
        }
        tempobs.push_back(total / (double)(conf.xo()));
        _corrfunc.push_back(tempobs);
      } 
    }
    else
    {
      // average over configuration
      std::vector<double> temp(max-initial, 0.0);
      for(int p=0; p<(conf.xo()); p++)
      {
        int refspin = conf[initial];
        int index = 0;
        for(int i=initial; i<max; i++)
        {
          temp[index] += refspin*conf[i];
          index++;
        }
        conf.propagate();
      }
      // push value for correlation function estimator
      for(unsigned int i=0; i< temp.size(); i++)
        (_corrfunc[i]).push_back(temp[i]/(double)(conf.xo())); 
    }
  }
 
  void PROBE::output_to_file(const std::string& ofile)
  {
    std::ofstream output;
    output.open(ofile);
    output << std::setw(104) << std::setfill('-') << "-" << '\n';
    output << std::setfill(' ');
    output << std::setw(99) << std::left << "system size:" << std::setw(5)  
           << _ns << '\n';
    output << std::setw(99) << std::left << "d-field:"     << std::setw(5)
           << _df << '\n';
    output << std::setw(99) << std::left << "beta:"        << std::setw(5)
           << _bt << '\n';
    output << std::setw(95) << std::left << "boundary conditions:";
    if(_bc) output << std::setw(9) << "periodic" << '\n';
    else    output << std::setw(9) << "open"     << '\n';
    output << std::setw(104) << std::setfill('-') << "-" << '\n';
    output << std::setfill(' ');
    output << std::setw(78) << std::left << "energy:" 
           << std::setw(18) << std::setprecision(3) 
           << std::fixed << _energy.ave() 
           << std::setw(8) << std::setprecision(1) 
           << std::scientific << _energy.err()  
           << '\n';
    output << std::setw(78) << std::left << "specific heat:" 
           << std::setw(18) << std::setprecision(3) 
           << std::fixed 
           << (_nsquared.ave() - _nave.ave()*_nave.ave() - _nave.ave())/6 
           << std::setw(8) << std::setprecision(1) 
           << std::scientific << _energy.err()  
           << '\n';
    output << std::setw(78) << std::left << "staggered mag. susceptibility:" 
           << std::setw(18) << std::setprecision(5) 
           << std::fixed << _smsusc.ave()/_bt
           << std::setw(8) << std::setprecision(1) 
           << std::scientific << _smsusc.err()/_bt  
           << '\n';
    output << std::setw(78) << std::left << "staggered mag. fluctuations:" 
           << std::setw(18) << std::setprecision(5) 
           << std::fixed << _smfluc.ave() 
           << std::setw(8) << std::setprecision(1) 
           << std::scientific << _smfluc.err()  
           << '\n';
    output << std::setw(78) << std::left << "string order susceptibility:" 
           << std::setw(18) << std::setprecision(5) 
           << std::fixed << _sosusc.ave()/_bt
           << std::setw(8) << std::setprecision(1) 
           << std::scientific << _sosusc.err()/_bt  
           << '\n';
    output << std::setw(78) << std::left << "string order fluctuations:" 
           << std::setw(18) << std::setprecision(5) 
           << std::fixed << _sofluc.ave() 
           << std::setw(8) << std::setprecision(1) 
           << std::scientific << _sofluc.err()  
           << '\n';
    output << std::setw(78) << std::left << "staggered magnetization:" 
           << std::setw(18) << std::setprecision(5) 
           << std::fixed << _sm.ave()
           << std::setw(8) << std::setprecision(1) 
           << std::scientific << _sm.err()
           << '\n';
    output << std::setw(78) << std::left << "string order:" 
           << std::setw(18) << std::setprecision(5) 
           << std::fixed << _so.ave() 
           << std::setw(8) << std::setprecision(1) 
           << std::scientific << _so.err()  
           << '\n';
    output << std::setw(104) << std::setfill('-') << "-" << '\n';
    output << std::setfill(' ');
    output << std::setw(6)  << "Site";
    if(!_corrfunc.empty()){
      output << std::setw(14) << "<S_i^zS_j^z>"
             << std::setw(10) << "error";
    }
    output << '\n';
    output << std::setw(104) << std::setfill('-') << "-" << '\n';
    output << std::setfill(' ');
    for(unsigned int i=0; i<_corrfunc.size(); i++)
    {
      output << std::setw(6)  << i;
      if(!_corrfunc.empty()){
        if(i<_corrfunc.size()){
          output << std::setw(14) << std::fixed  
                                  << std::setprecision(3) 
                                  << _corrfunc[i].ave()
                 << std::setw(10) << std::scientific
                                  << std::setprecision(1)
                                  << _corrfunc[i].err();
        }
        else{
          output << std::setw(15) << " " 
                 << std::setw(10) << " ";
        }
      }
      output << '\n';
    }                      
  }         
}
