// Written by Peter Kutz.

#ifndef CONSTANTS_H
#define CONSTANTS_H

#ifndef __MINMAX_DEFINED
#  define max(a,b)    (((a) > (b)) ? (a) : (b))
#  define min(a,b)    (((a) < (b)) ? (a) : (b))
#endif

#define LERP(a,b,t) (1-t)*a + t*b



// Don't modify the values of these here.
// Modify the values of these in Constants.cpp instead.
extern const int theMillisecondsPerFrame;
extern const int theDim[3];
extern const double theCellSize;
extern const double theAirDensity;
extern const double theBuoyancyAlpha;
extern const double theBuoyancyBeta;	
extern const double theBuoyancyAmbientTemperature;
extern const double theVorticityEpsilon;





#endif