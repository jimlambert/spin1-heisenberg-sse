#include <iostream>
#include "CONFIG.h"

using namespace std;

int main()
{
  int N = 10;
  SSE::CONFIG testconf(N, 0.1, 0.0, 0.1, true);

  //testconf.disp_wgts();
  //testconf.disp_extprbs();
  return 0;
}
