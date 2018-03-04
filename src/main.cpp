#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <boost/program_options.hpp>
#include "CONFIG.h"
#include "PROBE.h"
#include "OBSERVABLE.h"

namespace po = boost::program_options;
using namespace std;

int main(int ac, char *av[])
{
  int n, equil, simul;
  double t, d, e;
  bool pbc;
  string ofile;
  po::options_description desc("Available options");
  desc.add_options()
    ("default",                    "default program settings")
    ("help",                       "display options available")
    ("size",  po::value<int>(),    "system size")
    ("equil", po::value<int>(),    "number of equilibration steps")
    ("simul", po::value<int>(),    "number of simulation steps")
    ("temp",  po::value<double>(), "system temperature")
    ("dfld",  po::value<double>(), "control uniaxial anisotropy")
    ("eps",   po::value<double>(), "controls energy offset")
    ("pbc",   po::value<bool>(),   "controls periodic boundary conditions")
    ("of",    po::value<string>(), "output file (including directory)");
  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);

  n     = 4;
  equil = 10000;
  simul = 10000;
  t     = 0.1;
  d     = 0.0;
  e     = 0.0;
  pbc   = true;
  ofile = "defoutput.dat";
  if(vm.count("default"))
  {
    cout << "proceeding with default settings" << endl;
  }
  else
  {
    if(vm.count("help"))
    {
      cout << desc << endl;
      return 1;
    }
    if(vm.count("size")) n = vm["size"].as<int>();
    else
    {
      cout << "missing size, using default: " << n << endl;
    }  
    if(vm.count("equil")) equil = vm["equil"].as<int>();
    else
    {
      cout << "missing equilibration steps, using default: " << equil << endl;
    }
    if(vm.count("simul")) simul = vm["simul"].as<int>();
    else
    {
      cout << "missing simulation steps, using default: " << simul << endl;
    }
    if(vm.count("temp")) t = vm["temp"].as<double>();
    else
    {
      cout << "missing temperature, using default: " << t << endl;
    }
    if(vm.count("dfld")) d = vm["dfld"].as<double>();
    else
    {
      cout << "missing D-field, using default: " << d << endl;
    }
    if(vm.count("eps")) e = vm["eps"].as<double>();
    else
    {
      cout << "missing epsilon, using default: " << e << endl;
    }
    if(vm.count("pbc")) pbc = vm["pbc"].as<bool>();
    else
    {
      cout << "missing boundary conditions, using default: " << pbc << endl;
    }
    if(vm.count("of")) ofile = vm["of"].as<string>();
    else
    {
      cout << "missing output file, using default: " << ofile << endl;
    }
  }
  
  SSE::CONFIG testconf(n, t, d, e, pbc);
  SSE::PROBE  probe(testconf, 
                    M_PI,
                    true,     // measure energy
                    true,     // measure square of string correlator
                    true,     // measure average of square of sum
                    true,     // measure S_z correlator
                    true,     // measure self correlations
                    true,     // measure spin averages
                    true      // measure string correlations
                   );
  //testconf.diagupdt();
  //testconf.disp_config();
  //testconf.loopupdt();
  //testconf.disp_config();
  //testconf.diagupdt();
  //testconf.loopupdt();
  //testconf.disp_config();
  
  for(int i=0; i<equil; i++)
  {
    testconf.diagupdt();
    testconf.loopupdt();
    testconf.expoupdt();
    cout << (int)((double)i/((double)equil)*100) << " %" << '\r' << flush;
  }
  cout << "equilibration complete" << endl;

  for(int i=0; i<simul; i++)
  {
    testconf.diagupdt();
    testconf.loopupdt();
    probe(testconf);
    cout << (int)((double)i/((double)simul)*100) << " %" << '\r' << flush;
  }
  cout << "simulation complete" << endl;
  probe.output_to_file(ofile);
  return 0;
}
