#include "CONFIG.h"

// VRTTABLE   -   allowed vertices for the spin-1 AFM Heisenberg model with
//                uniaxial anisotropy
// TYPTABLE   -   Contains array that indicates whether the vertex at the 
//                corresponding index is diagonal or off-diagonal
// OUTTABLE   -   contains output vertices each corresponding to exiting on a
//                given leg of the bare vertex
// PRBTABLE   -   contains the bounds for the transition probabilities given an
//                input leg, input vertex, and change type.

typedef static std::array<<std::array<int, 4>,    17>   VRTTABLE;   
typedef static std::array<bool,                   17>   TYPTABLE;
typedef static std::array<<std::array<int, 4>,    136>  OUTTABLE;   
typedef static std::array<<std::array<double, 4>, 136>  PRBTABLE;   


VRTTABLE verts = 
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
  { 0,  -1, -1,  0},  // 17
};

TYPETABLE types = 
{
  true;               // 1
  true;               // 2
  true;               // 3
  true;               // 4
  true;               // 5
  true;               // 6
  true;               // 7
  true;               // 8
  true;               // 9
  false;              // 10 
  false;              // 11
  false;              // 12
  false;              // 13
  false;              // 14
  false;              // 15
  false;              // 16
  false;              // 17
};  

OUTTABLE outvrts;
PRBTABLE extprbs;

namespace SSE
{
   
}
