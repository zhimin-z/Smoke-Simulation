#include "smoke_sim.h"
#include "constants.h"
#include "open_gl_headers.h"
#include "stb_image_write.h"
#include "custom_output.h"
#include "basic_math.h"

SmokeSim::SmokeSim() : mFrameNum(0), mTotalFrameNum(0), mRecordEnabled(false)
{
   reset();
}

SmokeSim::~SmokeSim()
{
}

void SmokeSim::reset()
{
   mGrid.reset();
	mTotalFrameNum = 0;
}

void SmokeSim::step()
{
	double dt = 0.4;
	mTotalFrameNum++;
  // TODO : Write the simulation step

	// Step0: Gather user forces
	mGrid.updateSources();

	// Step1: Calculate new velocities
	mGrid.advectVelocity(dt);
	mGrid.addExternalForces(dt);
	mGrid.project(dt);

	// Step2: Calculate new density 
	mGrid.advectDensity(dt);

	// Step3: Calculate new temperature
	mGrid.advectTemperature(dt);
}

void SmokeSim::setRecording(bool on, int width, int height)
{
   if (on && ! mRecordEnabled)  // reset counter
   {
      mFrameNum = 0;
   }
   mRecordEnabled = on;
	
	recordWidth = width;
	recordHeight = height;
}

bool SmokeSim::isRecording()
{
   return mRecordEnabled;
}

void SmokeSim::draw(const Camera& c)
{
   drawAxes(); 
   mGrid.draw(c);
   if (mRecordEnabled) grabScreen();
}

void SmokeSim::drawAxes()
{
	glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
		glDisable(GL_LIGHTING);

		glLineWidth(2.0); 
		glBegin(GL_LINES);
			glColor3f(1.0, 0.0, 0.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(1.0, 0.0, 0.0);

			glColor3f(0.0, 1.0, 0.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, 1.0, 0.0);

			glColor3f(0.0, 0.0, 1.0);
			glVertex3f(0.0, 0.0, 0.0);
			glVertex3f(0.0, 0.0, 1.0);
		glEnd();
	glPopAttrib();
}

void SmokeSim::grabScreen() // OVERHAULED BY PETER KUTZ. USING STB_IMAGE_WRITE INSTEAD OF DEVIL.
{
	
	
	if (mFrameNum > 9999) exit(0);
	
	// Save the smoke densities:
	char smoke_filename[2048];
#if WIN32
	sprintf_s(smoke_filename, 2048, "smoke_%04d.txt", mFrameNum); // Use snprintf if your compiler supports C99.
#elif __APPLE__
  sprintf(smoke_filename, "smoke_%04d.txt", mFrameNum);
#endif
	mGrid.saveSmoke(smoke_filename);

	
	// Save an image:
	unsigned char* bitmapData = new unsigned char[3 * recordWidth * recordHeight];

	for (int i=0; i<recordHeight; i++) 
	{
		glReadPixels(0,i,recordWidth,1,GL_RGB, GL_UNSIGNED_BYTE, 
			bitmapData + (recordWidth * 3 * ((recordHeight-1)-i)));
	}

	char anim_filename[2048];
#if WIN32
	sprintf_s(anim_filename, 2048, "smoke_%04d.png", mFrameNum); 
#elif __APPLE__
  sprintf(anim_filename, "smoke_%04d.png", mFrameNum);
#endif
	
	stbi_write_png(anim_filename, recordWidth, recordHeight, 3, bitmapData, recordWidth * 3);
	
	delete [] bitmapData;
	
	mFrameNum++;
}

int SmokeSim::getTotalFrames() {
	return mTotalFrameNum;
}
