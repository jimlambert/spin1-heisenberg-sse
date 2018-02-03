#include <iostream>
#include "CONFIG.h"

using namespace std;

int main()
{
  int N = 2;
  SSE::CONFIG testconf(N, 0.01, 0.0, 0.01, true);
  /*  
  testconf.diagupdt();
  testconf.loopupdt();
  testconf.disp_config();
  cout << "---------------------------------" << endl;
  testconf.diagupdt();
  testconf.loopupdt();
  testconf.disp_config();
  cout << "---------------------------------" << endl;
  testconf.diagupdt();
  testconf.loopupdt();
  testconf.disp_config();
  */
  double enrgy = 0.0;

  for(int i=0; i<20000; i++)
  {
    testconf.diagupdt();
    testconf.loopupdt();
    testconf.expoupdt();
    cout << (double)i/(double)20000 << '\r' << flush;
  }
  cout << endl;

  for(int i=0; i<20000; i++)
  {
    testconf.diagupdt();
    testconf.loopupdt();
    testconf.expoupdt();
    enrgy += testconf.no();
    cout << (double)i/(double)20000 << '\r' << flush;
  }

  cout << endl;

  cout << endl;
  cout << "energy: " << (-0.01*(enrgy/(20000*N)) + testconf.eo()) 
       << endl; 
  return 0;
}
