#include <iostream>
#include "CONFIG.h"

using namespace std;

int main()
{
  int N = 10;
  SSE::CONFIG testconf(N, 0.1, 0.0, 0.0, true);

  /*
  for(int i=0; i<N; i++)
  {
    cout << testconf[i] << '\t'; 
  }
  cout << endl;
  */
  return 0;
}
