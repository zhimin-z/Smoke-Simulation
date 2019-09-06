// ==========================================================================
// Copyright (C) 2009 Aline Normoyle
// ==========================================================================
#ifndef smokeSim_H_
#define smokeSim_H_

#include "mac_grid.h"

class Camera;
class SmokeSim
{
public:
   SmokeSim();
   virtual ~SmokeSim();

   virtual void reset();
   virtual void step();
   virtual void draw(const Camera& c);
   //virtual void setGridDimensions(int x, int y, int z); // REMOVED BY PETER KUTZ.
   virtual void setRecording(bool on, int width, int height);
   virtual bool isRecording();
	
	// PETER KUTZ:
	int getTotalFrames();

protected:
   virtual void drawAxes();
   virtual void grabScreen();

protected:
	MACGrid mGrid;
	bool mRecordEnabled;
	int mFrameNum;
	int mTotalFrameNum; // PETER KUTZ.
	
	// PETER KUTZ
	int recordWidth;
	int recordHeight;
};

#endif