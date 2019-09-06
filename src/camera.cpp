#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "open_gl_headers.h" // PETER KUTZ.
#include "basic_math.h" // PETER KUTZ.

#include "camera.h"

glm::dvec3 Camera::dfltEye = glm::dvec3(5.0, 15.0, -25.0);
glm::dvec3 Camera::dfltUp = glm::dvec3(0.0, 1.0, 0.0);
glm::dvec3 Camera::dfltLook = glm::dvec3(5.0, 10.0, 0.0);
float Camera::dfltVfov = 60.0;
float Camera::dfltAspect = 1.0;
float Camera::dfltNear = 0.5;
float Camera::dfltFar = 120.0; // Increased by PETER KUTZ.
float Camera::dfltSpeed = 0.1;
float Camera::dfltTurnRate = 1.0 * (BasicMath::PI/180.0);

Camera::Camera() 
{   
   myDir = NONE; myTurnDir = NONE;
   reset();
}

Camera::~Camera() {}

void Camera::reset()
{
   mSpeed = dfltSpeed;
   mTurnRate = dfltTurnRate;
   mVfov = dfltVfov;
   mAspect = dfltAspect;
   mNear = dfltNear;
   mFar = dfltFar;

   // Calculate the initial heading & pitch
   // Note that  eye[0] = radius*cos(h)*cos(p); and  eye[1] = radius*sin(p);
   mPitch = -std::asin(dfltEye[1]/glm::length(dfltEye));
   mHeading = std::acos(dfltEye[0]/(glm::length(dfltEye)*std::cos(mPitch)));
   //printf("INIT: %f %f\n", mPitch, mHeading);

   set(dfltEye, dfltLook, dfltUp);
}

void Camera::draw()
{
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glm::mat4 persp = glm::perspective(mVfov, mAspect, mNear, mFar);
   glLoadMatrixf(glm::value_ptr(persp));

   float m[16];
   m[0] = v[0]; m[4] = v[1]; m[8] = v[2];  m[12] = -glm::dot(eye, v); 
   m[1] = u[0]; m[5] = u[1]; m[9] = u[2];  m[13] = -glm::dot(eye, u); 
   m[2] = n[0]; m[6] = n[1]; m[10] = n[2]; m[14] = -glm::dot(eye, n); 
   m[3] = 0.0;  m[7] = 0.0;  m[11] = 0.0;  m[15] = 1.0;
   glMatrixMode(GL_MODELVIEW);
   glLoadMatrixf(m); 

   glGetDoublev(GL_MODELVIEW_MATRIX, glm::value_ptr(myModelMatrix));
   glGetDoublev(GL_PROJECTION_MATRIX, glm::value_ptr(myProjMatrix));
   glGetIntegerv(GL_VIEWPORT, glm::value_ptr(myViewport));
}

const glm::dvec3& Camera::getUp() const
{
   return u;
}

const glm::dvec3& Camera::getBackward() const
{
   return n;
}

const glm::dvec3& Camera::getRight() const
{
   return v;
}

glm::vec3 Camera::getRelativePosition(float left, float up, float forward)
{
   glm::vec3 direction = up * glm::vec3(u) + left * glm::vec3(v) - forward * glm::vec3(n);
   return glm::vec3(eye) + direction;  // Move along forward axis 
}

void Camera::getViewport(int& x, int& y, int& w, int& h)
{
   x = myViewport[0];
   y = myViewport[1];
   w = myViewport[2];
   h = myViewport[3];
}

void Camera::getProjection(
   float* vfov, float* aspect, float* zNear, float* zFar)
{
   *vfov = mVfov; *aspect = mAspect; *zNear = mNear; *zFar = mFar;
}

void Camera::setPosition(const glm::vec3& pos)
{
   eye = pos;
}

const glm::dvec3& Camera::getPosition() const
{
   return eye;
}

void Camera::setProjection(
   float vfov, float aspect, float zNear, float zFar)
{
   mVfov = vfov;
   mAspect = aspect;
   mNear = zNear;
   mFar = zFar;
}

float Camera::heading() const
{
   return mHeading;
}

float Camera::pitch() const
{
   return mPitch;
}

void Camera::set(const glm::dvec3& eyepos, const glm::dvec3& look, const glm::dvec3& up)
{
	eye = eyepos;
	n = eyepos - look;
	v = glm::cross(up,n);
	u = glm::cross(n,v);
	mRadius = glm::length(n); // cache this distance

	u = glm::normalize(u);
	v = glm::normalize(v);
	n = glm::normalize(n);
}

void Camera::print()
{
   std::cout << "EYE: " << eye.x << ", " << eye.y << ", " << eye.z << std::endl;
   std::cout << "RIGHT: " << v.x << ", " << v.y << ", " << v.z << std::endl;
   std::cout << "UP: " << u.x << ", " << u.y << ", " << u.z << std::endl;
   std::cout << "N: " << n.x << ", " << n.y << ", " << n.z << std::endl;

   printf("-----------------------\n");
}

void Camera::move(float dV, float dU, float dN)
{
   eye += glm::dvec3(dU * glm::vec3(u) + dV * glm::vec3(v) + dN * glm::vec3(n));
}

void Camera::orbit(float h, float p)
{
  //printf("PITCH: %f\n", p);
  //printf("HEADING: %f\n", h);
  //printf("RADIUS: %f\n", mRadius);

   glm::dvec3 rotatePt; // Calculate new location around sphere having mRadius
   rotatePt[0] = double(mRadius) * std::cos(h)*std::cos(p);
   rotatePt[1] = double(mRadius) * std::sin(p);
   rotatePt[2] = double(mRadius) * std::sin(h)*std::cos(p);

   glm::dvec3 lookAt = eye - n * double(mRadius);
   set(lookAt - rotatePt, lookAt /* look */, glm::dvec3(0.0, 1.0, 0.0) /* up Approx */);
}

void Camera::orbitLeft(int scale) 
{
   myTurnDir = TL;
   mHeading -= mTurnRate * scale; // Inverted by PETER KUTZ.
   orbit(mHeading, pitch());
}

void Camera::moveLeft(int scale) // => move along v
{    
   myDir = L;
   move(-mSpeed*scale, 0.0, 0.0);
}

void Camera::orbitRight(int scale)
{
   myTurnDir = TR;
   mHeading += mTurnRate * scale; // Inverted by PETER KUTZ.
   orbit(mHeading, pitch());
}

void Camera::moveRight(int scale) // => move along v
{
   myDir = R;
   move(mSpeed*scale, 0.0, 0.0);   
}

void Camera::orbitUp(int scale)
{
   myTurnDir = TU; 
   mPitch = min(BasicMath::PI/2.0 - 0.01, mPitch + mTurnRate * scale);
   orbit(heading(), mPitch);
}

void Camera::moveUp(int scale) // => move along +u
{
   myDir = U;
   move(0.0, mSpeed*scale, 0.0);   
}

void Camera::orbitDown(int scale)
{
   myTurnDir = TD; 
   mPitch = max(-BasicMath::PI/2.0 + 0.01, mPitch - mTurnRate*scale);
   orbit(heading(), mPitch);
}

void Camera::moveDown(int scale) // => move along -u
{
   myDir = D;
   move(0.0, -mSpeed*scale, 0.0);   
}

void Camera::moveForward(int scale) // => move along -n
{
   myDir = F; 
   move(0.0, 0.0, -mSpeed*scale);      
   mRadius += -mSpeed*scale;  // Also "zoom" into radius
}

void Camera::moveBack(int scale) // => move along n
{
   myDir = B; 
   move(0.0, 0.0, mSpeed*scale);   
   mRadius += mSpeed*scale;  // Also "zoom" out radius
}

void Camera::turn(glm::dvec3& v1, glm::dvec3& v2, float amount)
{
   double cosTheta = std::cos(amount);
   double sinTheta = std::sin(amount);

   float vX =  cosTheta*v1[0] + sinTheta*v2[0]; 
   float vY =  cosTheta*v1[1] + sinTheta*v2[1]; 
   float vZ =  cosTheta*v1[2] + sinTheta*v2[2]; 

   float nX = -sinTheta*v1[0] + cosTheta*v2[0]; 
   float nY = -sinTheta*v1[1] + cosTheta*v2[1]; 
   float nZ = -sinTheta*v1[2] + cosTheta*v2[2]; 

   v1.x = vX, v1.y = vY, v1.z = vZ;
   v2.x = nX, v2.y = nY, v2.z = nZ;
}

void Camera::turnLeft(int scale) // rotate around u
{
   myTurnDir = TL; 
   turn(v, n, -mTurnRate*scale);
}

void Camera::turnRight(int scale) // rotate around u
{
   myTurnDir = TR;
   turn(v, n, mTurnRate*scale);
}

void Camera::turnUp(int scale) // rotate around v
{
   myTurnDir = TU; 
   turn(n, u, mTurnRate*scale);
}

void Camera::turnDown(int scale) // rotate around v
{
   myTurnDir = TD; 
   turn(n, u, -mTurnRate*scale);
}

bool Camera::screenToWorld(int screenX, int screenY, glm::vec3& worldCoords)
{
   double x, y, z;
   worldCoords = glm::unProject(glm::vec3(screenX, screenY, 0.0), glm::mat4(myModelMatrix), glm::mat4(myProjMatrix), glm::vec4(myViewport));
   return true;
}

bool Camera::worldToScreen(const glm::vec3& worldCoords, int& screenX, int& screenY)
{
   double x, y, z;
   glm::vec3 screenCoords = glm::project(worldCoords, glm::mat4(myModelMatrix), glm::mat4(myProjMatrix), glm::vec4(myViewport));

   screenX = (int) screenCoords.x;
   screenY = (int) screenCoords.y;
   return true;
}

glm::mat4 Camera::cameraToWorldMatrix()
{
   glm::mat4 tmp;
   tmp = myModelMatrix;
   tmp = glm::inverse(tmp);
   return tmp;
}
