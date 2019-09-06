// Written by Peter Kutz.

#include "constants.h"

const int theMillisecondsPerFrame = 10;

#ifdef _DEBUG
const int theDim[3] = {10, 10, 5};
#else
const int theDim[3] = {40,40,20};
#endif

const double theCellSize = 1.0;
const double theAirDensity = 1.0;

const double theBuoyancyAlpha =0; // Gravity's effect on the smoke particles.
const double theBuoyancyBeta = 0; // Buoyancy's effect due to temperature difference.	
const double theBuoyancyAmbientTemperature = 0; // Ambient temperature.

const double theVorticityEpsilon = 2.0;