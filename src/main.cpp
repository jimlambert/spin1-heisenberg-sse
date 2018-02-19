#include <iostream>
#include <cmath>
#include "CONFIG.h"
#include "PROBE.h"
#include "OBSERVABLE.h"

using namespace std;

int main()
{
  int N = 2;
  SSE::CONFIG testconf(N, 1, 0.0, 0.01, true);
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
  testconf.diagupdt();
  testconf.loopupdt();
  testconf.disp_config(); 
  std::cout << " --- " << std::endl;
  testconf.diagupdt();
  testconf.loopupdt();
  testconf.disp_config(); 
  std::cout << " --- " << std::endl;
  /*
  for(int i=0; i<10000; i++)
  {
    testconf.diagupdt();
    testconf.loopupdt();
    testconf.expoupdt();
    cout << (double)i/((double)10000) << '\r' << flush;
  }
  cout << "equilibration complete" << endl;

  for(int i=0; i<10000; i++)
  {
    testconf.diagupdt();
    testconf.loopupdt();
    testconf.expoupdt();
    probe(testconf);
    cout << (double)i/((double)10000) << '\r' << flush;
  }
  cout << "simulation complete" << endl;
  probe.output_to_file("./testoutput.dat");
  */
  return 0;
}
