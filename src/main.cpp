// CIS 563 Smoke Simulation : Semi-Lagrangian Smoke Simulator
// Updated by Harmony Li, Copyright (c) 2015 University of Pennsylvania

#include "main.h"

//--------------------------
//-------- MAIN ------------
//--------------------------

int main(int argc, char** argv) {
  // Initialize GLFW Window and OpenGL
  if (init(argc, argv)) {
    isRunning = true;

    // GLFW main loop
    mainLoop();
  }

  glfwTerminate();
  exit(1);
}

bool init(int argc, char** argv) {
  glfwSetErrorCallback(errorCallback);
  if (!glfwInit()) {
    return false;
  }

  width = 800;
  height = 800;

  savedWidth = width;
  savedHeight = height;

  window = glfwCreateWindow(width, height, "CIS 563 Smoke Simulation", NULL, NULL);
  if (!window) {
    glfwTerminate();
    return false;
  }

  glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, keyCallback);
  glfwSetCursorPosCallback(window, cursorCallback);
  glfwSetWindowSizeCallback(window, windowResizeCallback);

  // Set up GL context
  glewExperimental = GL_TRUE;
  if (glewInit() != GLEW_OK) {
    return false;
  }

  initCamera();
  setGL();
  
  return true;
}

void initCamera()
{
   double w = theDim[0]*theCellSize;   
   double h = theDim[1]*theCellSize;   
   double d = theDim[2]*theCellSize;   
   double angle = 0.5*theCamera.dfltVfov*BasicMath::PI/180.0;
   double dist;
   if (w > h) dist = w*0.5/std::tan(angle);  // aspect is 1, so i can do this
   else dist = h*0.5/std::tan(angle);
   theCamera.dfltEye.x = w*0.5, theCamera.dfltEye.y = h*0.5, theCamera.dfltEye.z = -(dist+d*0.5);
   theCamera.dfltLook.x = w*0.5, theCamera.dfltLook.y = h*0.5, theCamera.dfltLook.z = 0.0;
   theCamera.reset();
}

void setGL() {
  glClearColor(0.1, 0.1, 0.1, 1.0);

  glEnable(GL_BLEND);
  glEnable(GL_ALPHA_TEST);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);
  glShadeModel(GL_SMOOTH);

  glEnable(GL_NORMALIZE);
  glDisable(GL_LIGHTING);
  glCullFace(GL_BACK);
}

void mainLoop() {
  while(!glfwWindowShouldClose(window)) {
    glfwPollEvents();

    if (isRunning) {
      theSmokeSim.step();
    }

    // Keep track of time
	  theFpsTracker.timestamp();

    // Draw elements
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	  theCamera.draw();
	  theSmokeSim.draw(theCamera);
	  drawOverlay();
    glfwSwapBuffers(window);
  }

  glfwDestroyWindow(window);
  glfwTerminate();
}

void drawOverlay()
{
  // Draw Overlay
  glColor4f(1.0, 1.0, 1.0, 1.0);
  glPushAttrib(GL_LIGHTING_BIT);
  glDisable(GL_LIGHTING);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glRasterPos2f(0.01, 0.01);
     
  char info[1024];
  sprintf(info, "CIS 563 Smoke Simulation | Framerate: %3.1f  |  Frame: %u  |  %s", 
    theFpsTracker.fpsAverage(), theSmokeSim.getTotalFrames(),
    theSmokeSim.isRecording()? "Recording..." : "");
 
  glfwSetWindowTitle(window, info);

  glPopAttrib();
}

// CLEAN UP
void shutDown(int returnCode) {
#ifdef __APPLE__
  glfwTerminate();
#endif
  exit(returnCode);
}

void errorCallback(int error, const char* message) {
  fputs(message, stderr);
}

void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
  if (action == GLFW_PRESS) {
    if (mods & GLFW_MOD_SHIFT) {
      switch(key) {
      case GLFW_KEY_PERIOD:
        isRunning = true;
        break;
      case GLFW_KEY_COMMA:
        theSmokeSim.reset();
        break;
      }
    } else {
      switch (key) {
      case GLFW_KEY_0:
        MACGrid::theRenderMode = MACGrid::CUBES;
        break;
      case GLFW_KEY_1:
        MACGrid::theRenderMode = MACGrid::SHEETS;
        break;
      case GLFW_KEY_V:
        MACGrid::theDisplayVel = !MACGrid::theDisplayVel;
        break;
      case GLFW_KEY_SPACE:
        theCamera.reset();
        break;
      case GLFW_KEY_R:
        theSmokeSim.setRecording(!theSmokeSim.isRecording(), savedWidth, savedHeight);
        break;
      case GLFW_KEY_EQUAL:
        isRunning = false;
        break;      
      case GLFW_KEY_ESCAPE:
        glfwSetWindowShouldClose(window, GL_TRUE);
        break;
      }
    }
  }
}

void cursorCallback(GLFWwindow* window, double x, double y) {
  int deltaX = (int)(lastX - x);
  int deltaY = (int)(lastY - y);
  bool moveLeftRight = abs(deltaX) > abs(deltaY);
  bool moveUpDown = !moveLeftRight;

  if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
    if (moveLeftRight && deltaX > 0) theCamera.orbitLeft(deltaX);
    else if (moveLeftRight && deltaX < 0) theCamera.orbitRight(-deltaX);
    else if (moveUpDown && deltaY > 0) theCamera.orbitUp(deltaY);
    else if (moveUpDown && deltaY < 0) theCamera.orbitDown(-deltaY);
  } else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS) {
    if (moveUpDown && deltaY > 0) theCamera.moveForward(deltaY);
    else if (moveUpDown && deltaY < 0) theCamera.moveBack(-deltaY);
  }

  if (glfwGetKey(window, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS || glfwGetKey(window, GLFW_KEY_LEFT_ALT) == GLFW_PRESS) {
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
      if (moveLeftRight && deltaX > 0) theCamera.moveLeft(deltaX);
      else if (moveLeftRight && deltaX < 0) theCamera.moveRight(-deltaX);
      else if (moveUpDown && deltaY > 0) theCamera.moveUp(deltaY);
      else if (moveUpDown && deltaY < 0) theCamera.moveDown(-deltaY);
    }
  }

  lastX = x;
  lastY = y;
}

void windowResizeCallback(GLFWwindow* window, int w, int h) {
  width = w, savedWidth = w;
  height = h, savedHeight = h;

  // Update viewport
  glViewport(0, 0, width, height);

  // Update camera projection's aspect ratio
  float vfov, aspect, zNear, zFar;
  theCamera.getProjection(&vfov, &aspect, &zNear, &zFar);
  theCamera.setProjection(vfov, ((GLfloat) width)/height, zNear, zFar);
}
