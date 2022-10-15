# Smoke+Fire Simulation

Implemented a 3D Semi-Lagrangian smoke fluid simulation using the Navier-Stokes equations and a marker-and-cell (MAC) grid in C++ and rendered in real-time using OpenGL for [CIS 563 - Phys Based Animation](https://www.coursicle.com/penn/courses/CIS/563/) 2015 Spring Project 4. Velocity, density, and temperature are advected onto the smoke particles and buoyancy and vorticity confinement external forces are applied on the system. Incompressibility is maintained by doing a pressure solve using a conjugate gradient solver and then projecting it onto the system.

See https://youtu.be/LQhsG7FdlMU
