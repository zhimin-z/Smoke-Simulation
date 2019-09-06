// ==========================================================================
// Copyright (C) 2008 Aline Normoyle
// ==========================================================================

#ifndef camera_H_
#define camera_H_

#include "open_gl_headers.h" // PETER KUTZ.

#include <algorithm>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <iostream>

#include "constants.h"

class Camera
{
public:
   Camera();
   virtual ~Camera();

   // Draw projection and eyepoint
   virtual void draw();

   // Print eyepoint position and basis
   virtual void print();

   // Initialize the camera with glyLookAt parameters
   virtual void set(const glm::dvec3& eyepos, const glm::dvec3& look, const glm::dvec3& up);

   // Get camera state
   virtual void setPosition(const glm::vec3& pos);
   virtual const glm::dvec3& getPosition() const;
   virtual const glm::dvec3& getUp() const;
   virtual const glm::dvec3& getBackward() const;
   virtual const glm::dvec3& getRight() const;
   virtual glm::vec3 getRelativePosition(float left, float up, float forward);
   virtual float heading() const;
   virtual float pitch() const;

   // Camera frustrum managements
   virtual void setProjection(
      float vfov, float aspect, float zNear, float zFar);
   virtual void getProjection(
      float* vfov, float* aspect, float* zNear, float* zFar);
   virtual void getViewport(int& x, int& y, int& w, int& h);

   // Relative movement commands
   virtual void moveLeft(int scale = 1.0);
   virtual void moveRight(int scale = 1.0);
   virtual void moveUp(int scale = 1.0);
   virtual void moveDown(int scale = 1.0);
   virtual void moveForward(int scale = 1.0);
   virtual void moveBack(int scale = 1.0);

   virtual void turnLeft(int scale = 1.0);
   virtual void turnRight(int scale = 1.0);
   virtual void turnUp(int scale = 1.0);
   virtual void turnDown(int scale = 1.0);

   virtual void orbitLeft(int scale = 1.0);
   virtual void orbitRight(int scale = 1.0);
   virtual void orbitUp(int scale = 1.0);
   virtual void orbitDown(int scale = 1.0);

   // Reset to original state
   virtual void reset();

   // Conversion utilities between screen and world coordinates
   virtual bool screenToWorld(int screenX, int screenY, glm::vec3& worldCoords);
   virtual bool worldToScreen(const glm::vec3& worldCoords, int& screenX, int& screenY);

   // Get camera to world matrix
   virtual glm::mat4 cameraToWorldMatrix();

protected:
   enum Dir { NONE, F, B, L, R, U, D, TL, TR, TU, TD} myDir, myTurnDir;
   virtual void turn(glm::dvec3& v, glm::dvec3& n, float amount);
   virtual void move(float dU, float dV, float dN);
   virtual void orbit(float h, float p);

protected:
   float mSpeed, mTurnRate;

   glm::dvec3 eye; // camera position
   float mHeading, mPitch, mRadius;
   float mVfov, mAspect, mNear, mFar; // projection parameters
   
   // Basis of camera local coord system
   glm::dvec3 u; // up
   glm::dvec3 v; // v points right
   glm::dvec3 n; // -n points forward

   // Cache useful values
   glm::dmat4 myProjMatrix;
   glm::dmat4 myModelMatrix;
   glm::ivec4 myViewport;

public:

   // Defaults
   static glm::dvec3 dfltEye, dfltUp, dfltLook;
   static float dfltVfov, dfltAspect, dfltNear, dfltFar; 
   static float dfltSpeed, dfltTurnRate;
};

#endif
