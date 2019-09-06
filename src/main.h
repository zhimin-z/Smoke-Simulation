// CIS 563 Smoke Simulation : Semi-Lagrangian Smoke Simulator
// Updated by Harmony Li, Copyright (c) 2015 University of Pennsylvania

#ifndef MAIN_H
#define MAIN_H

#include "open_gl_headers.h"

#include <stdio.h>
#include <cmath>
#include <glm/glm.hpp>

#include "basic_math.h"
#include "camera.h"
#include "constants.h"
#include "fps.h"
#include "smoke_sim.h"

// Geometry and whatnot
SmokeSim theSmokeSim;
Camera theCamera;
mmc::FpsTracker theFpsTracker;

// Window parameters
GLFWwindow *window;
int width;
int height;
int savedWidth = 0;
int savedHeight = 0;

// UI Helpers
double lastX = 0, lastY = 0;
int theButtonState = 0;
int theModifierState = 0;
bool isRunning = true;

//-------------------------
//------ SETUP STUFF ------
//-------------------------
bool init(int, char**);
void initCamera();
void setGL();

//------------------------
//---- GLFW CALLBACKS ----
//------------------------
void mainLoop();
void drawOverlay();

void errorCallback(int, const char*);
void keyCallback(GLFWwindow*, int, int, int, int);
void cursorCallback(GLFWwindow*, double, double);
void windowResizeCallback(GLFWwindow*, int, int);


#endif /* MAIN_H */