// ==========================================================================
// Copyright (C) 2008 Aline Normoyle
// ==========================================================================

#ifndef GridData_H_
#define GridData_H_

#pragma warning(disable: 4244 4267 4996)
#include <iostream>
#include <algorithm>
#include <glm/glm.hpp>
#include <vector> // PETER KUTZ.
#include "Constants.h"// PETER KUTZ.

#define SIGN(x) (x > 0 ? +1 : (x < 0 ? -1 : 0))

// GridData is capable of storing any data in a grid
// Columns are indexed with i and increase with increasing x
// Rows are indexed with j and increase with z
// Stacks are indexed with k and increase with y
//
// GridData is initialized by global variables theDim and theCellSize
// defined in GridData.cpp.  theDim defines the number of cells in 
// each X,Y, Z direction.  theCellSize defines the size of each cell.
// GridData's world space dimensions extend from (0,0,0) to mMax, where mMax is
// (theCellSize*theDim[0], theCellSize*theDim[1], theCellSize*theDim[2])
class GridData
{
public:
   GridData();
   GridData(const GridData& orig);
   virtual ~GridData();
   virtual GridData& operator=(const GridData& orig);

   // Initialize underlying data structure with dlftValue
   virtual void initialize(double dfltValue = 0.0);

   // Returns editable data at index (i,j,k).
   // E.g. to set data on this object, call mygriddata(i,j,k) = newval
   virtual double& operator()(int i, int j, int k);
   virtual const double operator()(int i, int j, int k) const; // PETER KUTZ.

   // Given a point in world coordinates, return the corresponding
   // value from this grid. mDfltValue is returned for points
   // outside of our grid dimensions
   virtual double interpolate(const glm::dvec3& pt);
   double GridData::HCInterpolate(double q_1, double q0, double q1, double q2, double x);
   // Access underlying data structure (for use with other UBLAS objects)
   std::vector<double>& data();

   // Given a point in world coordinates, return the cell index (i,j,k)
   // corresponding to it
   virtual void getCell(const glm::dvec3& pt, int& i, int& j, int& k);

protected:

   virtual glm::dvec3 worldToSelf(const glm::dvec3& pt) const;
   std::vector<double> mData;//
   double mDfltValue;//
   glm::dvec3 mMax;//

};

class GridDataX : public GridData
{
public:
   GridDataX();
   virtual ~GridDataX();
   virtual void initialize(double dfltValue = 0.0);
   virtual double& operator()(int i, int j, int k);
   virtual const double operator()(int i, int j, int k) const;
   virtual glm::dvec3 worldToSelf(const glm::dvec3& pt) const;
};

class GridDataY : public GridData
{
public:
   GridDataY();
   virtual ~GridDataY();
   virtual void initialize(double dfltValue = 0.0);
   virtual double& operator()(int i, int j, int k);
   virtual const double operator()(int i, int j, int k) const;
   virtual glm::dvec3 worldToSelf(const glm::dvec3& pt) const;
};

class GridDataZ : public GridData
{
public:
   GridDataZ();
   virtual ~GridDataZ();
   virtual void initialize(double dfltValue = 0.0);
   virtual double& operator()(int i, int j, int k);
   virtual const double operator()(int i, int j, int k) const;
   virtual glm::dvec3 worldToSelf(const glm::dvec3& pt) const;
};

#endif
