/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
//#define GHOST
#define REACTIVE
#pragma once
#include "core.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <deque>
#include <queue>
#include <ctime>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <signal.h>
#include "boost/thread.hpp"
#include "Matrix.hpp"
#include "MatrixOperations.hpp"

#define USE_GL

#define MATERIALS 2
#define DIMENSIONS 2
#define GC 2

#include "grid.hpp"

template<Processor P>
struct Mesh {
  typedef Mesh2D<real, P, NUMBER_VARIABLES> type;
};

template<Processor P>
struct LevelSet {
  typedef Mesh2D<real, P, 2> type;
};

StridedCell<real, NUMBER_VARIABLES> typedef Cell;

/*///////////////////////////////////////////////////////////////////////////////////////////
// Read input flow field variables
int read() {
double var[1];
std::ifstream inputfile;
inputfile.open("INPUT.txt");

while (! inputfile.eof()){
    inputfile >> var[0];
}
inputfile.close();
return 0;
};
///////////////////////////////////////////////////////////////////////////////////////////*/


const real X1 = 10.0;
const real X2 = 30.0;
const real PCJ = 21.54;
const real rho_CJ = 1.795;
const real pShock = 2.0 * PCJ;
const real rhoShock = 1.0 * rho_CJ;
const real pAir = 1.0;
const real rhoAir = 1.0;
const real Tol = 1.0e-6;
const real m1 = 0.0;
const real cell = 5;

// Sinusoidal perturbation
const real PI = 3.14159265359;
const real AMP = 5.0;           // Amplitude
const real LAMBDA = 10.0;       // Wavelength


// Put in a 90 degree opening
//const real XCORNER = 1250.0;
//const real YCORNER = 240.0;
//const real GA = 0.04;
//const real SPACING = 40.0;


/*
// Put in a rectnagular block as a small perturbation
const real block_x1 = 1400.0;
const real block_x2 = 1450.0;
const real block_y1 = 280.0;
const real block_y2 = 285.0;
*/

//rough tube

const real number = 20; //number of the cubes
const real length_x = 15; //horizontal length of the cubes
const real length_y = 15; //vertical length of the cubes
const real space = 60; //space between the cubes in a row
const real phase = 30; //phase shift between top and bot cubes
const real start_x = 300; //the start point of the first cube
const real diameter = 300; //diameter of the tube

/*
const real X11 = 300.0;
const real Y11 = 0.0;
const real X12 = 315.0;
const real Y12 = 15.0;
const real X21 = 330.0;
const real Y21 = 285.0;
const real X22 = 345.0;
const real Y22 = 300.0;
const real X31 = 360.0;
const real Y31 = 0.0;
const real X32 = 375.0;
const real Y32 = 15.0;
const real X41 = 390.0;
const real Y41 = 285.0;
const real X42 = 405.0;
const real Y42 = 300.0;
const real X51 = 420.0;
const real Y51 = 0.0;
const real X52 = 435.0;
const real Y52 = 15.0;
const real X61 = 450.0;
const real Y61 = 285.0;
const real X62 = 465.0;
const real Y62 = 300.0;
const real X71 = 480.0;
const real Y71 = 0.0;
const real X72 = 495.0;
const real Y72 = 15.0;
const real X81 = 505.0;
const real Y81 = 285.0;
const real X82 = 520.0;
const real Y82 = 300.0;
const real X91 = 535.0;
const real Y91 = 0.0;
const real X92 = 550.0;
const real Y92 = 15.0;
const real X101 = 565.0;
const real Y101 = 285.0;
const real X102 = 580.0;
const real Y102 = 300.0;
const real X103 = 595.0;
const real Y103 = 0.0;
const real X104 = 610.0;
const real Y104 = 15.0;
const real X105 = 625.0;
const real Y105 = 285;
const real X106 = 640.0;
const real Y106 = 300.0;
const real X107 = 655.0;
const real Y107 = 0;
const real X108 = 670.0;
const real Y108 = 15.0;
const real X109 = 685.0;
const real Y109 = 285;
const real X110 = 700.0;
const real Y110 = 300.0;
const real X111 = 715.0;
const real Y111 = 0;
const real X112 = 730.0;
const real Y112 = 15.0;
const real X113= 745.0;
const real Y113 = 285;
const real X114 = 760.0;
const real Y114 = 300.0;
const real X115 = 775.0;
const real Y115 = 0;
const real X116 = 790.0;
const real Y116 = 15.0;
const real X117 = 805.0;
const real Y117 = 285;
const real X118 = 820.0;
const real Y118 = 300.0;
const real X119 = 835.0;
const real Y119 = 0;
const real X120 = 850.0;
const real Y120 = 15.0;
const real X121 = 865.0;
const real Y121 = 285;
const real X122 = 880.0;
const real Y122 = 300.0;
*/

// Parameters for 2-step kinetics (unstable H2/O2)
const real Q = 21.365;
const real TS = 5.0373;
const real EI = 5.414 * TS;
const real ER = 1.0 * TS;
const real KI = 1.0022;
const real KR = 4.0;
const real specific_heat_ratio = 1.32;

/*// Parameters for 2-step kinetics (unstable C2H2+2.5O2, p0=100kpa)
const real Q = 75.77;
const real TS = 7.65;
const real EI = 5.4 * TS;
const real ER = 1.0 * TS;
const real KI = 0.8316;
const real KR = 4;
const real specific_heat_ratio = 1.19;
*/

/*// Parameters for 2-step kinetics (stable C2H2/O2/Ar)
const real Q = 19.7;
const real TS = 7.6051;
const real EI = 4.8 * TS;
const real ER = 1.0 * TS;
const real KI = 0.139;
const real KR = 0.2;
const real specific_heat_ratio = 1.212;
*/

// Parameters for single-step Arrhenius kinetics (unstable H2/O2)
//const real Q = 50;
//const real KR = 85.5; // 
//const real Ea = 30.0;
//const real KR = 16.45;
//const real KR = 7.25;
//const real KR = 3.55;
//const real Ea = 20.0;
//const real KR = 8.5;
//const real KR = 1.6;
//const real Ea = 15.0;
//const real KR = 80.2;
//const real Ea = 30.0;
//const real specific_heat_ratio = 1.2;

/*real Counter = 2.0;
real Check = 1.0;
const real Plot_step = 0.5;
const int Skip_lines = 5;
real Plot_counter = 1.0;*/

real Counter = 1.0;
real Check = 1.0;
//const int len_inert = 500; //in terms of cells (multiplied by resolution)
real density_new = 1.0;
const real Delta_rho = 0.0;
const real frame_interface = 1.0;
const real Plot_step = 9.5;
real Plot_counter = 1.0;
const int Skip_lines = 5;

__device__ __host__ __forceinline__ real gamma(const Cell u) {
  return specific_heat_ratio;
}
__device__ __host__ __forceinline__ real gamma() {
  return specific_heat_ratio;
}

__device__ __host__ __forceinline__ real p0(const Cell u) {
  return 0.0;
  //return 0.87e8;
}
__device__ __host__ __forceinline__ real p0() {
  return 0.0;
  //return 0.87e8;
}

#include "boundaryconditions.hpp"
#include "flux.hpp"
#include "wavespeed.hpp"
#include "HLLC.hpp"
#include "Solver.hpp"
#include "initialconditions.hpp"
#include "render.hpp"
#include "opengl.hpp"
#ifdef GHOST
#include "ghostfluid.hpp"
#endif
#ifdef REACTIVE
#include "source.hpp"
#include "shockdetect.hpp"
#endif

struct ImageOutputs {
  std::string prefix;
  int plotVariable;
  ColourMode colourMode;
  real min;
  real max;
};

#include "SDF/BoundaryMesh.hpp"
#include "SDF/Polyhedron.cu"
#include "SDF/ConnectedEdge.cu"
#include "SDF/Edge.cu"
#include "SDF/Face.cu"
#include "SDF/Vertex.cu"
#include "SDF/ConnectedFace.cu"
#include "SDF/ConnectedVertex.cu"
#include "SDF/ScanConvertiblePolygon.cu"
#include "SDF/ScanConvertiblePolyhedron.cu"
#include "SDF/BoundaryMesh.cu"

#include "kernels/boundaryconditions.ipp"
#include "kernels/flux.ipp"
#ifdef GHOST
#include "kernels/ghostfluid.ipp"
#endif
#ifdef REACTIVE
#include "kernels/source.ipp"
#include "kernels/shockdetect.ipp"
#endif
#include "kernels/HLLC.ipp"
#include "kernels/wavespeed.ipp"
