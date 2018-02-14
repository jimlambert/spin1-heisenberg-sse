#include <vector>
#include <complex>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "OBSERVABLE.h"
#include "CONFIG.h"
#include "PROBE.h"

namespace SSE
{
  void PROBE::meas_energy(CONFIG& conf)
  {
    double val = (-1.0/_bt)*(conf.no()) + _ns*_offset;
    _energy.push_back(val);
  }

  void PROBE::meas_sqarstrg(CONFIG& conf)
  {
    int max, initial;
    if(_bc){
      max = _ns/2;
      initial = 0;
    }
    else{
      max = _ns/2;
      initial = 0;
    } 
    double total = 0.0; 
    std::complex<double> i(0, 1);
    for(int p=0; p<(conf.xo()); p++)
    {
      for(int k=initial; k<max; k++)
      {
        int sum = 0.0;
        for(int j=k; j<max; j++) sum += conf[j];
        double x = (double)conf[k]*std::exp(i*M_PI*(double)sum).real()
                  *(double)conf[max-1];
        total += x*x;
      }
      conf.propagate();
    }
    _sqarstrg.push_back(total / (double)(conf.xo()));
  }
 
  void PROBE::meas_sqarspin(CONFIG& conf)
  {
    double total = 0.0;
    std::complex<double> i(0, 1);
    for(int p=0; p<(conf.xo()); p++)
    {
      for(int k=0; k<_ns; k++)
      {
        total += std::exp(i*_kval*(double)k).real()*conf[k];
      }
      conf.propagate();
    }
    total = total/(double)conf.xo();
    _sqarspin.push_back(total*total);
  }

  void PROBE::meas_corrfunc(CONFIG& conf)
  {
    int max;
    int initial;
    if(_bc){  // PBC
      initial = 0;
      max = _ns;
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

  void PROBE::meas_sqarcorr(CONFIG& conf)
  {
    int max;
    int initial;
    if(_bc){  // PBC
      initial = 0;
      max = _ns;
    } 
    else{     // OBC
      initial = _ns /2;
      max = _ns;
    }   
    if(_sqarcorr.empty())
    {
      for(int i=initial; i<max; i++)
      {
        OBSERVABLE<double> tempobs;
        double total = 0.0;
        for(int p=0; p<(conf.xo()); p++)
        {
          total += conf[i]*conf[i];
          conf.propagate();
        }
        tempobs.push_back(total /(double)(conf.xo()));
        _sqarcorr.push_back(tempobs);
      }
    }
    else
    {
      std::vector<double> temp(max-initial, 0.0);
      for(int p=0; p<(conf.xo()); p++)
      {
        int index = 0;
        for(int i=initial; i<max; i++)
        {
          temp[index] += conf[i] * conf[i];
          index++; 
        }
        conf.propagate();
      } 
      for(int i=initial; i<max; i++)
        _sqarcorr[i].push_back(temp[i]/(double)(conf.xo()));
    }
  }
 
  void PROBE::meas_spinavrg(CONFIG& conf)
  {
    int max;
    int initial;
    if(_bc){  // PBC
      initial = 0;
      max = _ns;
    } 
    else{     // OBC
      initial = _ns /2;
      max = _ns;
    }
    if(_spinavrg.empty())
    {
      std::vector<double>  temp(max-initial, 0);
      for(int p=0; p<conf.xo(); p++)
      {
        int index = 0;
        for(int i=initial; i<max; i++)
        {
          temp[index] += conf[i];
          index++;
        }
        conf.propagate();
      }
      for(unsigned int i=0; i<temp.size(); i++)
      {
        OBSERVABLE<double> tempobs;
        tempobs.push_back(temp[i]/(double)conf.xo());
        _spinavrg.push_back(tempobs);
      }
    }
    else
    {
      std::vector<double> temp(max-initial, 0.0);
      for(int p=0; p<conf.xo(); p++)
      {
        int index=0;
        for(int i=initial; i<max; i++)
        {
          temp[index] += conf[i];
          index++;
        }
        conf.propagate();
      }
      for(unsigned int i=0; i<temp.size(); i++)
        _spinavrg[i].push_back(temp[i]/(double)conf.xo());
    }     
  }

  void PROBE::meas_strgcorr(CONFIG& conf)
  {
    int max, initial;
    if(_bc){
      max = _ns / 2;
      initial = 0;
    }
    else{
      max = _ns;
      initial = _ns / 2;
    }   
    if(_strgcorr.empty())
    {
      std::vector<double> temp(max-initial, 0);
      std::complex<double> i(0, 1);
      for(int p=0; p<(conf.xo()); p++)
      {
        int index = 0;
        for(int k=initial; k<max; k++)
        {
          double sum = 0.0;
          for(int j=k; j<max; j++) sum += conf[j];
          temp[index] += (double)conf[k]*std::exp(i*M_PI*(double)sum).real()
                        *(double)conf[max];
          index++;
        }
        conf.propagate();
      }
      for(int i=0; i<max; i++)
      {
        OBSERVABLE<double> tempobs;
        tempobs.push_back(temp[i]/(double)(conf.xo()));
        _strgcorr.push_back(tempobs);
      }
    }
    else
    {
      std::complex<double> i(0, 1);
      for(int p=0; p<(conf.xo()); p++)
      {
        for(int k=initial; k<max; k++)
        {
          double sum = 0.0;
          for(int j=k; j<max; j++) sum += conf[j];
          _strgcorr[k].push_back((double)conf[k]*std::exp(i*M_PI*sum).real()
                                *(double)conf[max]); 
        } 
        conf.propagate();
      }
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
    output << std::setw(88) << std::left << "energy:" 
           << std::setw(8) << std::setprecision(3) 
           << std::fixed << _energy.ave() 
           << std::setw(8) << std::setprecision(1) 
           << std::scientific << _energy.err()  
           << '\n';
    output << std::setw(88) << std::left << "string corr. squared:"
           << std::setw(8) << std::setprecision(3) 
           << std::fixed << _sqarstrg.ave() 
           << std::setw(8) << std::setprecision(1) 
           << std::scientific << _sqarstrg.err()
           << '\n';
    output << std::setw(88) << "staggered mag. squared:" 
           << std::setw(8) << std::setprecision(3) 
           << std::fixed << _sqarspin.ave()
           << std::setw(8) << std::setprecision(1) 
           << std::scientific << _sqarspin.err()
           << '\n';
    output << std::setw(104) << std::setfill('-') << "-" << '\n';
    output << std::setfill(' ');
    output << std::setw(6)  << "Site";
    if(!_corrfunc.empty()){
      output << std::setw(14) << "<S_i^zS_j^z>"
             << std::setw(10) << "error";
    }
    if(!_sqarcorr.empty()){
      output << std::setw(13) << "<(S_i^z)^2>"
             << std::setw(10) << "error";
    }
    if(!_spinavrg.empty()){
      output << std::setw(9)  << "<S_i^z>"
             << std::setw(10) << "error";
    }
    if(!_strgcorr.empty()){
      output << std::setw(22) << "String Correlations"
             << std::setw(10) << "error";
    }
    output << '\n';
    output << std::setw(104) << std::setfill('-') << "-" << '\n';
    output << std::setfill(' ');
    for(unsigned int i=0; i<_spinavrg.size(); i++)
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
      if(!_sqarcorr.empty()){
        if(i<_sqarcorr.size()){
          output << std::setw(13) << std::fixed
                                  << std::setprecision(3)
                                  << _sqarcorr[i].ave()
                 << std::setw(10)  << std::scientific
                                  << std::setprecision(1)
                                  << _sqarcorr[i].err();
        }
        else{
          output << std::setw(13) << " " 
                 << std::setw(10)  << " ";
        }
      }
      if(!_spinavrg.empty()){
        if(i<_spinavrg.size()){
          output << std::setw(9) << std::fixed
                                 << std::setprecision(3)
                                 << _spinavrg[i].ave()
                 << std::setw(10) << std::scientific
                                 << std::setprecision(1)
                                 << _spinavrg[i].err();
        }
        else{
          output << std::setw(9) << " " 
                 << std::setw(10) << " ";
        }
      }
      if(!_strgcorr.empty()){
        if(i<_strgcorr.size()){
          output << std::setw(22) << std::fixed
                                  << std::setprecision(3)
                                  << _strgcorr[i].ave()
                 << std::setw(10) << std::scientific
                                  << std::setprecision(1)
                                  << _strgcorr[i].err();
        }
        else{
          output << std::setw(22) << " " 
                 << std::setw(10) << " ";
        }
      }
      output << '\n';
    }                      
  }         
}
