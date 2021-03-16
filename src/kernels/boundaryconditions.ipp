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
/*
  if (XDIR && u.exists(k, 0) && k < u.i(XCORNER)) {
    const int j = u.j(YCORNER);
    for (int n = 0; n < 2; n++) {
      u(k, j + n ) = u(k, j - n - 1);
      u(k, j + n , YMOMENTUM) = -u(k, j - n -1, YMOMENTUM);
    }
  }
  else if (!XDIR && u.exists(0, k) && k > u.j(YCORNER)) {
    const int i = u.i(XCORNER);
    for (int n = 0; n < 2; n++) {
      u(i - n - 1, k) = u(i + n, k);
      u(i - n - 1, k, XMOMENTUM) = -u(i + n, k, XMOMENTUM);
    }
  }*/
//////////
/*  if (XDIR && u.exists(k, 0) && k < u.i(XCORNER1)) {
    const int j = u.j(YCORNER1);
    for (int n = 0; n < 2; n++) {
      u(k, j - n - 1) = u(k, j + n);
      u(k, j - n - 1, YMOMENTUM) = -u(k, j + n, YMOMENTUM);

    }
  }
  else if (!XDIR && u.exists(0, k) && k < u.j(YCORNER1)) {
    const int i = u.i(XCORNER1);
    for (int n = 0; n < 2; n++) {
      u(i - n - 1, k) = u(i + n, k);
      u(i - n - 1, k, XMOMENTUM) = -u(i + n, k, XMOMENTUM);
    }
  }
*/
///////
  if (XDIR && u.exists(k, 0) && k > u.i(X11) && k < u.i(X12)) {

   /* const int j_y11 = u.j(Y11);
    for (int n = 0; n < 2; n++) {
      u(k, j_y11 + n + 1) = u(k, j_y11 - n);
      u(k, j_y11 + n + 1, YMOMENTUM) = -u(k, j_y11 + n + 1, YMOMENTUM);
    }*/

    const int j_y12 = u.j(Y12);
    for (int n = 0; n < 2; n++) {
      u(k, j_y12 - n - 1) = u(k, j_y12 + n);
      u(k, j_y12 - n - 1, YMOMENTUM) = -u(k, j_y12 - n - 1, YMOMENTUM);
    }

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y11) && k < u.j(Y12)) {

    const int i_x11 = u.i(X11);
    for (int n = 0; n < 2; n++) {
      u(i_x11 + n + 1, k) = u(i_x11 - n, k);
      u(i_x11 + n + 1, k, XMOMENTUM) = -u(i_x11 + n + 1, k, XMOMENTUM);
    }

    const int i_x12 = u.i(X12);
    for (int n = 0; n < 2; n++) {
      u(i_x12 - n - 1, k) = u(i_x12 + n, k);
      u(i_x12 - n - 1, k, XMOMENTUM) = -u(i_x12 - n - 1, k, XMOMENTUM);
    }

  }
 ////
   if (XDIR && u.exists(k, 0) && k > u.i(X21) && k < u.i(X22)) {

    const int j_y21 = u.j(Y21);
    for (int n = 0; n < 2; n++) {
      u(k, j_y21 + n + 1) = u(k, j_y21 - n);
      u(k, j_y21 + n + 1, YMOMENTUM) = -u(k, j_y21 + n + 1, YMOMENTUM);
    }

    /*const int j_y22 = u.j(Y22);
    for (int n = 0; n < 2; n++) {
      u(k, j_y22 - n - 1) = u(k, j_y22 + n);
      u(k, j_y22 - n - 1, YMOMENTUM) = -u(k, j_y22 - n - 1, YMOMENTUM);
    }*/

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y21) && k < u.j(Y22)) {

    const int i_x21 = u.i(X21);
    for (int n = 0; n < 2; n++) {
      u(i_x21 + n + 1, k) = u(i_x21 - n, k);
      u(i_x21 + n + 1, k, XMOMENTUM) = -u(i_x21 + n + 1, k, XMOMENTUM);
    }

    const int i_x22 = u.i(X22);
    for (int n = 0; n < 2; n++) {
      u(i_x22 - n - 1, k) = u(i_x22 + n, k);
      u(i_x22 - n - 1, k, XMOMENTUM) = -u(i_x22 - n - 1, k, XMOMENTUM);
    }

  }
  
  /////
    if (XDIR && u.exists(k, 0) && k > u.i(X31) && k < u.i(X32)) {

    /*const int j_y31 = u.j(Y31);
    for (int n = 0; n < 2; n++) {
      u(k, j_y31 + n + 1) = u(k, j_y31 - n);
      u(k, j_y31 + n + 1, YMOMENTUM) = -u(k, j_y31 + n + 1, YMOMENTUM);
    }*/

    const int j_y32 = u.j(Y32);
    for (int n = 0; n < 2; n++) {
      u(k, j_y32 - n - 1) = u(k, j_y32 + n);
      u(k, j_y32 - n - 1, YMOMENTUM) = -u(k, j_y32 - n - 1, YMOMENTUM);
    }

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y31) && k < u.j(Y32)) {

    const int i_x31 = u.i(X31);
    for (int n = 0; n < 2; n++) {
      u(i_x31 + n + 1, k) = u(i_x31 - n, k);
      u(i_x31 + n + 1, k, XMOMENTUM) = -u(i_x31 + n + 1, k, XMOMENTUM);
    }

    const int i_x32 = u.i(X32);
    for (int n = 0; n < 2; n++) {
      u(i_x32 - n - 1, k) = u(i_x32 + n, k);
      u(i_x32 - n - 1, k, XMOMENTUM) = -u(i_x32 - n - 1, k, XMOMENTUM);
    }

  }
  
  ////
    if (XDIR && u.exists(k, 0) && k > u.i(X41) && k < u.i(X42)) {

    const int j_y41 = u.j(Y41);
    for (int n = 0; n < 2; n++) {
      u(k, j_y41 + n + 1) = u(k, j_y41 - n);
      u(k, j_y41 + n + 1, YMOMENTUM) = -u(k, j_y41 + n + 1, YMOMENTUM);
    }

   /* const int j_y42 = u.j(Y42);
    for (int n = 0; n < 2; n++) {
      u(k, j_y42 - n - 1) = u(k, j_y42 + n);
      u(k, j_y42 - n - 1, YMOMENTUM) = -u(k, j_y42 - n - 1, YMOMENTUM);
    }*/

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y41) && k < u.j(Y42)) {

    const int i_x41 = u.i(X41);
    for (int n = 0; n < 2; n++) {
      u(i_x41 + n + 1, k) = u(i_x41 - n, k);
      u(i_x41 + n + 1, k, XMOMENTUM) = -u(i_x41 + n + 1, k, XMOMENTUM);
    }

    const int i_x42 = u.i(X42);
    for (int n = 0; n < 2; n++) {
      u(i_x42 - n - 1, k) = u(i_x42 + n, k);
      u(i_x42 - n - 1, k, XMOMENTUM) = -u(i_x42 - n - 1, k, XMOMENTUM);
    }

  }
  
  ////
    if (XDIR && u.exists(k, 0) && k > u.i(X51) && k < u.i(X52)) {

   /* const int j_y51 = u.j(Y51);
    for (int n = 0; n < 2; n++) {
      u(k, j_y51 + n + 1) = u(k, j_y51 - n);
      u(k, j_y51 + n + 1, YMOMENTUM) = -u(k, j_y51 + n + 1, YMOMENTUM);
    }*/

    const int j_y52 = u.j(Y52);
    for (int n = 0; n < 2; n++) {
      u(k, j_y52 - n - 1) = u(k, j_y52 + n);
      u(k, j_y52 - n - 1, YMOMENTUM) = -u(k, j_y52 - n - 1, YMOMENTUM);
    }

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y51) && k < u.j(Y52)) {

    const int i_x51 = u.i(X51);
    for (int n = 0; n < 2; n++) {
      u(i_x51 + n + 1, k) = u(i_x51 - n, k);
      u(i_x51 + n + 1, k, XMOMENTUM) = -u(i_x51 + n + 1, k, XMOMENTUM);
    }

    const int i_x52 = u.i(X52);
    for (int n = 0; n < 2; n++) {
      u(i_x52 - n - 1, k) = u(i_x52 + n, k);
      u(i_x52 - n - 1, k, XMOMENTUM) = -u(i_x52 - n - 1, k, XMOMENTUM);
    }

  }
  ///
  if (XDIR && u.exists(k, 0) && k > u.i(X61) && k < u.i(X62)) {

    const int j_y61 = u.j(Y61);
    for (int n = 0; n < 2; n++) {
      u(k, j_y61 + n + 1) = u(k, j_y61 - n);
      u(k, j_y61 + n + 1, YMOMENTUM) = -u(k, j_y61 + n + 1, YMOMENTUM);
    }

    /*const int j_y62 = u.j(Y62);
    for (int n = 0; n < 2; n++) {
      u(k, j_y62 - n - 1) = u(k, j_y62 + n);
      u(k, j_y62 - n - 1, YMOMENTUM) = -u(k, j_y62 - n - 1, YMOMENTUM);
    }*/

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y61) && k < u.j(Y62)) {

    const int i_x61 = u.i(X61);
    for (int n = 0; n < 2; n++) {
      u(i_x61 + n + 1, k) = u(i_x61 - n, k);
      u(i_x61 + n + 1, k, XMOMENTUM) = -u(i_x61 + n + 1, k, XMOMENTUM);
    }

    const int i_x62 = u.i(X62);
    for (int n = 0; n < 2; n++) {
      u(i_x62 - n - 1, k) = u(i_x62 + n, k);
      u(i_x62 - n - 1, k, XMOMENTUM) = -u(i_x62 - n - 1, k, XMOMENTUM);
    }

  }
  ////
  if (XDIR && u.exists(k, 0) && k > u.i(X71) && k < u.i(X72)) {

   /* const int j_y71 = u.j(Y71);
    for (int n = 0; n < 2; n++) {
      u(k, j_y71 + n + 1) = u(k, j_y71 - n);
      u(k, j_y71 + n + 1, YMOMENTUM) = -u(k, j_y71 + n + 1, YMOMENTUM);
    }*/

    const int j_y72 = u.j(Y72);
    for (int n = 0; n < 2; n++) {
      u(k, j_y72 - n - 1) = u(k, j_y72 + n);
      u(k, j_y72 - n - 1, YMOMENTUM) = -u(k, j_y72 - n - 1, YMOMENTUM);
    }

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y71) && k < u.j(Y72)) {

    const int i_x71 = u.i(X71);
    for (int n = 0; n < 2; n++) {
      u(i_x71 + n + 1, k) = u(i_x71 - n, k);
      u(i_x71 + n + 1, k, XMOMENTUM) = -u(i_x71 + n + 1, k, XMOMENTUM);
    }

    const int i_x72 = u.i(X72);
    for (int n = 0; n < 2; n++) {
      u(i_x72 - n - 1, k) = u(i_x72 + n, k);
      u(i_x72 - n - 1, k, XMOMENTUM) = -u(i_x72 - n - 1, k, XMOMENTUM);
    }

  }
  ///
  if (XDIR && u.exists(k, 0) && k > u.i(X81) && k < u.i(X82)) {

    const int j_y81 = u.j(Y81);
    for (int n = 0; n < 2; n++) {
      u(k, j_y81 + n + 1) = u(k, j_y81 - n);
      u(k, j_y81 + n + 1, YMOMENTUM) = -u(k, j_y81 + n + 1, YMOMENTUM);
    }

    /*const int j_y82 = u.j(Y82);
    for (int n = 0; n < 2; n++) {
      u(k, j_y82 - n - 1) = u(k, j_y82 + n);
      u(k, j_y82 - n - 1, YMOMENTUM) = -u(k, j_y82 - n - 1, YMOMENTUM);
    }*/

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y81) && k < u.j(Y82)) {

    const int i_x81 = u.i(X81);
    for (int n = 0; n < 2; n++) {
      u(i_x81 + n + 1, k) = u(i_x81 - n, k);
      u(i_x81 + n + 1, k, XMOMENTUM) = -u(i_x81 + n + 1, k, XMOMENTUM);
    }

    const int i_x82 = u.i(X82);
    for (int n = 0; n < 2; n++) {
      u(i_x82 - n - 1, k) = u(i_x82 + n, k);
      u(i_x82 - n - 1, k, XMOMENTUM) = -u(i_x82 - n - 1, k, XMOMENTUM);
    }

  }
  //
  if (XDIR && u.exists(k, 0) && k > u.i(X91) && k < u.i(X92)) {

   /* const int j_y91 = u.j(Y91);
    for (int n = 0; n < 2; n++) {
      u(k, j_y91 + n + 1) = u(k, j_y91 - n);
      u(k, j_y91 + n + 1, YMOMENTUM) = -u(k, j_y91 + n + 1, YMOMENTUM);
    }*/

    const int j_y92 = u.j(Y92);
    for (int n = 0; n < 2; n++) {
      u(k, j_y92 - n - 1) = u(k, j_y92 + n);
      u(k, j_y92 - n - 1, YMOMENTUM) = -u(k, j_y92 - n - 1, YMOMENTUM);
    }

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y91) && k < u.j(Y92)) {

    const int i_x91 = u.i(X91);
    for (int n = 0; n < 2; n++) {
      u(i_x91 + n + 1, k) = u(i_x91 - n, k);
      u(i_x91 + n + 1, k, XMOMENTUM) = -u(i_x91 + n + 1, k, XMOMENTUM);
    }

    const int i_x92 = u.i(X92);
    for (int n = 0; n < 2; n++) {
      u(i_x92 - n - 1, k) = u(i_x92 + n, k);
      u(i_x92 - n - 1, k, XMOMENTUM) = -u(i_x92 - n - 1, k, XMOMENTUM);
    }

  }
  /////////////
  if (XDIR && u.exists(k, 0) && k > u.i(X101) && k < u.i(X102)) {

    const int j_y101 = u.j(Y101);
    for (int n = 0; n < 2; n++) {
      u(k, j_y101 + n + 1) = u(k, j_y101 - n);
      u(k, j_y101 + n + 1, YMOMENTUM) = -u(k, j_y101 + n + 1, YMOMENTUM);
    }

    /*const int j_y102 = u.j(Y102);
    for (int n = 0; n < 2; n++) {
      u(k, j_y102 - n - 1) = u(k, j_y102 + n);
      u(k, j_y102 - n - 1, YMOMENTUM) = -u(k, j_y102 - n - 1, YMOMENTUM);
    }*/

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y101) && k < u.j(Y102)) {

    const int i_x101 = u.i(X101);
    for (int n = 0; n < 2; n++) {
      u(i_x101 + n + 1, k) = u(i_x101 - n, k);
      u(i_x101 + n + 1, k, XMOMENTUM) = -u(i_x101 + n + 1, k, XMOMENTUM);
    }

    const int i_x102 = u.i(X102);
    for (int n = 0; n < 2; n++) {
      u(i_x102 - n - 1, k) = u(i_x102 + n, k);
      u(i_x102 - n - 1, k, XMOMENTUM) = -u(i_x102 - n - 1, k, XMOMENTUM);
    }

  }
  ////////////////
    if (XDIR && u.exists(k, 0) && k > u.i(X103) && k < u.i(X104)) {

    const int j_y104 = u.j(Y104);
    for (int n = 0; n < 2; n++) {
      u(k, j_y104 + n + 1) = u(k, j_y104 - n);
      u(k, j_y104 + n + 1, YMOMENTUM) = -u(k, j_y104 + n + 1, YMOMENTUM);
    }

    /*const int j_y102 = u.j(Y102);
    for (int n = 0; n < 2; n++) {
      u(k, j_y102 - n - 1) = u(k, j_y102 + n);
      u(k, j_y102 - n - 1, YMOMENTUM) = -u(k, j_y102 - n - 1, YMOMENTUM);
    }*/

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y103) && k < u.j(Y104)) {

    const int i_x103 = u.i(X103);
    for (int n = 0; n < 2; n++) {
      u(i_x103 + n + 1, k) = u(i_x103 - n, k);
      u(i_x103 + n + 1, k, XMOMENTUM) = -u(i_x103 + n + 1, k, XMOMENTUM);
    }

    const int i_x104 = u.i(X104);
    for (int n = 0; n < 2; n++) {
      u(i_x104 - n - 1, k) = u(i_x104 + n, k);
      u(i_x104 - n - 1, k, XMOMENTUM) = -u(i_x104 - n - 1, k, XMOMENTUM);
    }

  }
  ///////////////////////
    if (XDIR && u.exists(k, 0) && k > u.i(X105) && k < u.i(X106)) {

    const int j_y105 = u.j(Y105);
    for (int n = 0; n < 2; n++) {
      u(k, j_y105 + n + 1) = u(k, j_y105 - n);
      u(k, j_y105 + n + 1, YMOMENTUM) = -u(k, j_y105 + n + 1, YMOMENTUM);
    }

    /*const int j_y102 = u.j(Y102);
    for (int n = 0; n < 2; n++) {
      u(k, j_y102 - n - 1) = u(k, j_y102 + n);
      u(k, j_y102 - n - 1, YMOMENTUM) = -u(k, j_y102 - n - 1, YMOMENTUM);
    }*/

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y105) && k < u.j(Y106)) {

    const int i_x105 = u.i(X105);
    for (int n = 0; n < 2; n++) {
      u(i_x105 + n + 1, k) = u(i_x105 - n, k);
      u(i_x105 + n + 1, k, XMOMENTUM) = -u(i_x105 + n + 1, k, XMOMENTUM);
    }

    const int i_x106 = u.i(X106);
    for (int n = 0; n < 2; n++) {
      u(i_x106 - n - 1, k) = u(i_x106 + n, k);
      u(i_x106 - n - 1, k, XMOMENTUM) = -u(i_x106 - n - 1, k, XMOMENTUM);
    }

  }
  ///////////////////////
    if (XDIR && u.exists(k, 0) && k > u.i(X107) && k < u.i(X108)) {

    const int j_y108 = u.j(Y108);
    for (int n = 0; n < 2; n++) {
      u(k, j_y108 + n + 1) = u(k, j_y108 - n);
      u(k, j_y108 + n + 1, YMOMENTUM) = -u(k, j_y108 + n + 1, YMOMENTUM);
    }

    /*const int j_y102 = u.j(Y102);
    for (int n = 0; n < 2; n++) {
      u(k, j_y102 - n - 1) = u(k, j_y102 + n);
      u(k, j_y102 - n - 1, YMOMENTUM) = -u(k, j_y102 - n - 1, YMOMENTUM);
    }*/

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y107) && k < u.j(Y108)) {

    const int i_x107 = u.i(X107);
    for (int n = 0; n < 2; n++) {
      u(i_x107 + n + 1, k) = u(i_x107 - n, k);
      u(i_x107 + n + 1, k, XMOMENTUM) = -u(i_x107 + n + 1, k, XMOMENTUM);
    }

    const int i_x108 = u.i(X108);
    for (int n = 0; n < 2; n++) {
      u(i_x108 - n - 1, k) = u(i_x108 + n, k);
      u(i_x108 - n - 1, k, XMOMENTUM) = -u(i_x108 - n - 1, k, XMOMENTUM);
    }

  }
  ///////////////////////
    if (XDIR && u.exists(k, 0) && k > u.i(X109) && k < u.i(X110)) {

    const int j_y109 = u.j(Y109);
    for (int n = 0; n < 2; n++) {
      u(k, j_y109 + n + 1) = u(k, j_y109 - n);
      u(k, j_y109 + n + 1, YMOMENTUM) = -u(k, j_y109 + n + 1, YMOMENTUM);
    }

    /*const int j_y102 = u.j(Y102);
    for (int n = 0; n < 2; n++) {
      u(k, j_y102 - n - 1) = u(k, j_y102 + n);
      u(k, j_y102 - n - 1, YMOMENTUM) = -u(k, j_y102 - n - 1, YMOMENTUM);
    }*/

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y109) && k < u.j(Y110)) {

    const int i_x109 = u.i(X109);
    for (int n = 0; n < 2; n++) {
      u(i_x109 + n + 1, k) = u(i_x109 - n, k);
      u(i_x109 + n + 1, k, XMOMENTUM) = -u(i_x109 + n + 1, k, XMOMENTUM);
    }

    const int i_x110 = u.i(X110);
    for (int n = 0; n < 2; n++) {
      u(i_x110 - n - 1, k) = u(i_x110+ n, k);
      u(i_x110 - n - 1, k, XMOMENTUM) = -u(i_x110 - n - 1, k, XMOMENTUM);
    }

  }
  ///////////////////////
    if (XDIR && u.exists(k, 0) && k > u.i(X111) && k < u.i(X112)) {

    const int j_y112 = u.j(Y112);
    for (int n = 0; n < 2; n++) {
      u(k, j_y112 + n + 1) = u(k, j_y112 - n);
      u(k, j_y112 + n + 1, YMOMENTUM) = -u(k, j_y112+ n + 1, YMOMENTUM);
    }

    /*const int j_y102 = u.j(Y102);
    for (int n = 0; n < 2; n++) {
      u(k, j_y102 - n - 1) = u(k, j_y102 + n);
      u(k, j_y102 - n - 1, YMOMENTUM) = -u(k, j_y102 - n - 1, YMOMENTUM);
    }*/

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y111) && k < u.j(Y112)) {

    const int i_x111 = u.i(X111);
    for (int n = 0; n < 2; n++) {
      u(i_x111 + n + 1, k) = u(i_x111 - n, k);
      u(i_x111 + n + 1, k, XMOMENTUM) = -u(i_x111 + n + 1, k, XMOMENTUM);
    }

    const int i_x112 = u.i(X112);
    for (int n = 0; n < 2; n++) {
      u(i_x112 - n - 1, k) = u(i_x112+ n, k);
      u(i_x112 - n - 1, k, XMOMENTUM) = -u(i_x112 - n - 1, k, XMOMENTUM);
    }

  }
  ///////////////////////
    if (XDIR && u.exists(k, 0) && k > u.i(X113) && k < u.i(X114)) {

    const int j_y113 = u.j(Y113);
    for (int n = 0; n < 2; n++) {
      u(k, j_y113 + n + 1) = u(k, j_y113 - n);
      u(k, j_y113 + n + 1, YMOMENTUM) = -u(k, j_y113+ n + 1, YMOMENTUM);
    }

    /*const int j_y102 = u.j(Y102);
    for (int n = 0; n < 2; n++) {
      u(k, j_y102 - n - 1) = u(k, j_y102 + n);
      u(k, j_y102 - n - 1, YMOMENTUM) = -u(k, j_y102 - n - 1, YMOMENTUM);
    }*/

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y113) && k < u.j(Y114)) {

    const int i_x113 = u.i(X113);
    for (int n = 0; n < 2; n++) {
      u(i_x113 + n + 1, k) = u(i_x113 - n, k);
      u(i_x113+ n + 1, k, XMOMENTUM) = -u(i_x113+ n + 1, k, XMOMENTUM);
    }

    const int i_x114 = u.i(X114);
    for (int n = 0; n < 2; n++) {
      u(i_x114 - n - 1, k) = u(i_x114 + n, k);
      u(i_x114 - n - 1, k, XMOMENTUM) = -u(i_x114 - n - 1, k, XMOMENTUM);
    }

  }
   ///////////////////////
    if (XDIR && u.exists(k, 0) && k > u.i(X115) && k < u.i(X116)) {

    const int j_y116 = u.j(Y116);
    for (int n = 0; n < 2; n++) {
      u(k, j_y116 + n + 1) = u(k, j_y116 - n);
      u(k, j_y116 + n + 1, YMOMENTUM) = -u(k, j_y116+ n + 1, YMOMENTUM);
    }

    /*const int j_y102 = u.j(Y102);
    for (int n = 0; n < 2; n++) {
      u(k, j_y102 - n - 1) = u(k, j_y102 + n);
      u(k, j_y102 - n - 1, YMOMENTUM) = -u(k, j_y102 - n - 1, YMOMENTUM);
    }*/

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y115) && k < u.j(Y116)) {

    const int i_x115 = u.i(X115);
    for (int n = 0; n < 2; n++) {
      u(i_x115 + n + 1, k) = u(i_x115 - n, k);
      u(i_x115+ n + 1, k, XMOMENTUM) = -u(i_x115+ n + 1, k, XMOMENTUM);
    }

    const int i_x116 = u.i(X116);
    for (int n = 0; n < 2; n++) {
      u(i_x116 - n - 1, k) = u(i_x116 + n, k);
      u(i_x116 - n - 1, k, XMOMENTUM) = -u(i_x116 - n - 1, k, XMOMENTUM);
    }

  }
   ///////////////////////
    if (XDIR && u.exists(k, 0) && k > u.i(X117) && k < u.i(X118)) {

    const int j_y117 = u.j(Y117);
    for (int n = 0; n < 2; n++) {
      u(k, j_y117 + n + 1) = u(k, j_y117 - n);
      u(k, j_y117 + n + 1, YMOMENTUM) = -u(k, j_y117 + n + 1, YMOMENTUM);
    }

    /*const int j_y102 = u.j(Y102);
    for (int n = 0; n < 2; n++) {
      u(k, j_y102 - n - 1) = u(k, j_y102 + n);
      u(k, j_y102 - n - 1, YMOMENTUM) = -u(k, j_y102 - n - 1, YMOMENTUM);
    }*/

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y117) && k < u.j(Y118)) {

    const int i_x117 = u.i(X117);
    for (int n = 0; n < 2; n++) {
      u(i_x117 + n + 1, k) = u(i_x117 - n, k);
      u(i_x117 + n + 1, k, XMOMENTUM) = -u(i_x117 + n + 1, k, XMOMENTUM);
    }

    const int i_x118 = u.i(X118);
    for (int n = 0; n < 2; n++) {
      u(i_x118 - n - 1, k) = u(i_x118 + n, k);
      u(i_x118 - n - 1, k, XMOMENTUM) = -u(i_x118 - n - 1, k, XMOMENTUM);
    }

  }
   ///////////////////////
    if (XDIR && u.exists(k, 0) && k > u.i(X119) && k < u.i(X120)) {

    const int j_y120 = u.j(Y120);
    for (int n = 0; n < 2; n++) {
      u(k, j_y120 + n + 1) = u(k, j_y120 - n);
      u(k, j_y120 + n + 1, YMOMENTUM) = -u(k, j_y120+ n + 1, YMOMENTUM);
    }

    /*const int j_y102 = u.j(Y102);
    for (int n = 0; n < 2; n++) {
      u(k, j_y102 - n - 1) = u(k, j_y102 + n);
      u(k, j_y102 - n - 1, YMOMENTUM) = -u(k, j_y102 - n - 1, YMOMENTUM);
    }*/

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y119) && k < u.j(Y120)) {

    const int i_x119 = u.i(X119);
    for (int n = 0; n < 2; n++) {
      u(i_x119 + n + 1, k) = u(i_x119 - n, k);
      u(i_x119 + n + 1, k, XMOMENTUM) = -u(i_x119 + n + 1, k, XMOMENTUM);
    }

    const int i_x120 = u.i(X120);
    for (int n = 0; n < 2; n++) {
      u(i_x120 - n - 1, k) = u(i_x120 + n, k);
      u(i_x120 - n - 1, k, XMOMENTUM) = -u(i_x120 - n - 1, k, XMOMENTUM);
    }

  }
   ///////////////////////
    if (XDIR && u.exists(k, 0) && k > u.i(X121) && k < u.i(X122)) {

    const int j_y121 = u.j(Y121);
    for (int n = 0; n < 2; n++) {
      u(k, j_y121 + n + 1) = u(k, j_y121 - n);
      u(k, j_y121 + n + 1, YMOMENTUM) = -u(k, j_y121+ n + 1, YMOMENTUM);
    }

    /*const int j_y102 = u.j(Y102);
    for (int n = 0; n < 2; n++) {
      u(k, j_y102 - n - 1) = u(k, j_y102 + n);
      u(k, j_y102 - n - 1, YMOMENTUM) = -u(k, j_y102 - n - 1, YMOMENTUM);
    }*/

  }

  else if (!XDIR && u.exists(0, k) && k > u.j(Y121) && k < u.j(Y121)) {

    const int i_x121= u.i(X121);
    for (int n = 0; n < 2; n++) {
      u(i_x121 + n + 1, k) = u(i_x121 - n, k);
      u(i_x121 + n + 1, k, XMOMENTUM) = -u(i_x121 + n + 1, k, XMOMENTUM);
    }

    const int i_x122 = u.i(X122);
    for (int n = 0; n < 2; n++) {
      u(i_x122 - n - 1, k) = u(i_x122 + n, k);
      u(i_x122 - n - 1, k, XMOMENTUM) = -u(i_x122 - n - 1, k, XMOMENTUM);
    }

  }
}
  