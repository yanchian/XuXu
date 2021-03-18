/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#include "boundaryconditions.hpp"

template<BoundaryConditions BCs, bool XDIR, bool downstream>
__global__ void setBoundaryConditionsKernel(Mesh<GPU>::type u) {
  const bool YDIR = !XDIR, upstream = !downstream;

  const int i = blockIdx.x * blockDim.x + threadIdx.x - u.ghostCells();

  if (XDIR && u.exists(i, 0)) {
    for (int k = 0; k < NUMBER_VARIABLES; k++) {
      if (downstream) {
         u(i, u.activeNy() + 1, k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, u.activeNy() - 2, k);
        u(i, u.activeNy()    , k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, u.activeNy() - 1, k);
      // u(i, u.activeNy() + 1, k) = u(i, 1, k);
      // u(i, u.activeNy()    , k) =  u(i, 0, k);
      // u(i, u.activeNy() + 1, k) = u(i, u.activeNy() - 2, k); 
      // u(i, u.activeNy()    , k) = u(i, u.activeNy() - 1, k);
      //} else {
	  }
  }}
	  if (XDIR && u.exists(i, 0)) {
    for (int k = 0; k < NUMBER_VARIABLES; k++) {
      if (downstream) {
		  u(i, -2, k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, 1, k);
         u(i, -1, k) = (k == YMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(i, 0, k);
   //     u(i, -2, k) = u(i, u.activeNy() - 2, k); 
   //     u(i, -1, k) = u(i, u.activeNy() - 1, k);
       //u(i, -2, k) = u(i, 1, k);
       //u(i, -1, k) = u(i, 0, k);
      }
    }
  }
  if (!XDIR && u.exists(0, i)) {
    for (int k = 0; k < NUMBER_VARIABLES; k++) {
      if (downstream) {
         u(-2, i, k) = (k == XMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(1, i, k);
         u(-1, i, k) = (k == XMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(0, i, k);
        //u(-2, i, k) = u(1, i, k);
       // u(-1, i, k) = u(0, i, k); 
      } 
  }}
	  if (!XDIR && u.exists(0, i)) {
    for (int k = 0; k < NUMBER_VARIABLES; k++) {
      if (downstream) {
         u(u.activeNx() + 1, i, k) = (k == XMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(u.activeNx() - 2, i, k);
        u(u.activeNx()    , i, k) = (k == XMOMENTUM && BCs == REFLECTIVE ? -1.0 : 1.0) * u(u.activeNx() - 1, i, k);
        //u(u.activeNx() + 1, i, k) = u(u.activeNx() - 2, i, k);
       //u(u.activeNx()    , i, k) = u(u.activeNx() - 1, i, k);
      }
    }
  }
}

template<bool XDIR>
__global__ void setSpecialBoundaryConditionsKernel(Mesh<GPU>::type u) {
  const bool YDIR = !XDIR;

 // const int ki = blockIdx.x * blockDim.x + threadIdx.x - u.ghostCells();
 // const int kj = blockIdx.y * blockDim.y + threadIdx.y - u.ghostCells();
  const int k = blockIdx.x * blockDim.x + threadIdx.x - u.ghostCells();
  
for (int n_cube = 0; n_cube < number; n_cube++){
	
if (XDIR && u.exists(k, 0) && k > u.i(start_x + n_cube*space) && k < u.i(start_x + n_cube*space + length_x)) {

    const int j = u.j(length_y);
    for (int n = 0; n < 2; n++) {
      u(k, j - n - 1) = u(k, j + n);
      u(k, j - n - 1, YMOMENTUM) = -u(k, j - n - 1, YMOMENTUM);
    }

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(0) && k < u.j(length_y)) {

    const int i_1 = u.i(start_x + n_cube*space);
    for (int n = 0; n < 2; n++) {
      u(i_1 + n + 1, k) = u(i_1 - n, k);
      u(i_1 + n + 1, k, XMOMENTUM) = -u(i_1 + n + 1, k, XMOMENTUM);
    }

    const int i_2 = u.i(start_x + n_cube*space+length_x);
    for (int n = 0; n < 2; n++) {
      u(i_2 - n - 1, k) = u(i_2 + n, k);
      u(i_2 - n - 1, k, XMOMENTUM) = -u(i_2 - n - 1, k, XMOMENTUM);
    }
  }


else if (XDIR && u.exists(k, 0) && k > u.i(start_x + n_cube*space + phase) && k < u.i(start_x + n_cube*space + length_x + phase)) {

    const int j = u.j(start_x-length_y);
    for (int n = 0; n < 2; n++) {
      u(k, j + n + 1) = u(k, j - n);
      u(k, j + n + 1, YMOMENTUM) = -u(k, j + n + 1, YMOMENTUM);
    }

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(diameter - length_y) && k < u.j(diameter)) {

    const int i_1 = u.i(start_x + n_cube*space + phase);
    for (int n = 0; n < 2; n++) {
      u(i_1 + n + 1, k) = u(i_1 - n, k);
      u(i_1 + n + 1, k, XMOMENTUM) = -u(i_1 + n + 1, k, XMOMENTUM);
    }

    const int i_2 = u.i(start_x + n_cube*space+length_x + phase);
    for (int n = 0; n < 2; n++) {
      u(i_2 - n - 1, k) = u(i_2 + n, k);
      u(i_2 - n - 1, k, XMOMENTUM) = -u(i_2 - n - 1, k, XMOMENTUM);
    }
  }
}

}
