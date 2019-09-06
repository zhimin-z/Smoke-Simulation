#include "mac_grid.h"
#include "open_gl_headers.h"
#include "camera.h"
#include "custom_output.h"
#include "constants.h"

#include <math.h>
#include <map>
#include <stdio.h>

#undef max
#undef min
#include <fstream>
//#define __TEMP_FIELD__

// Globals
MACGrid target;


// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = false;

#define FOR_EACH_CELL \
	for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
	for(int j = 0; j < theDim[MACGrid::Y]; j++) \
	for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_CELL_REVERSE \
	for(int k = theDim[MACGrid::Z] - 1; k >= 0; k--)  \
	for(int j = theDim[MACGrid::Y] - 1; j >= 0; j--) \
	for(int i = theDim[MACGrid::X] - 1; i >= 0; i--) 

#define FOR_EACH_FACE \
	for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
	for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
	for(int i = 0; i < theDim[MACGrid::X]+1; i++) 

#define FOR_EACH_FACE_X \
	for (int k = 0; k < theDim[MACGrid::Z]; k++) \
	for (int j = 0; j < theDim[MACGrid::Y]; j++) \
	for (int i = 0; i < theDim[MACGrid::X]+1; i++)

#define FOR_EACH_FACE_Y \
	for (int k = 0; k < theDim[MACGrid::Z]; k++) \
	for (int j = 0; j < theDim[MACGrid::Y]+1; j++) \
	for (int i = 0; i < theDim[MACGrid::X]; i++)

#define FOR_EACH_FACE_Z \
	for (int k = 0; k < theDim[MACGrid::Z]+1; k++) \
	for (int j = 0; j < theDim[MACGrid::Y]; j++) \
	for (int i = 0; i < theDim[MACGrid::X]; i++)

MACGrid::MACGrid()
{
	initialize();
}

MACGrid::MACGrid(const MACGrid& orig)
{
	// TODO : Copy constructor for MAC Grid 
	mU = orig.mU;
	mV = orig.mV;
	mW = orig.mW;
	mP = orig.mP;
	mD = orig.mD;
	mT = orig.mT;
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
	// TODO : Copy constructor for MAC Grid
	if (&orig != this){
		mU = orig.mU;
		mV = orig.mV;
		mW = orig.mW;
		mP = orig.mP;
		mD = orig.mD;
		mT = orig.mT;   
	}
	return *this;
}

MACGrid::~MACGrid()
{
}

void MACGrid::reset()
{
	// TODO : Initialize the MAC Grid.
	mU.initialize(0);
	mV.initialize(0);
	mW.initialize(0);
	mP.initialize(0);
	mD.initialize(0);
	mT.initialize(0);
	calculateAMatrix();
}

void MACGrid::initialize()
{
	reset();
}

void MACGrid::updateSources()
{
	// TODO: Set initial values for density, temperature, velocity
	mV(2,1,1)=2.0;
	mV(1,2,1)=2.0;	
	mV(1,1,2)=2.0;
	mV(1,1,1)=2.0;

	mD(2,1,1)=1.0;
	mD(1,2,1)=1.0;	
	mD(1,1,2)=1.0;
	mD(1,1,1)=1.0;

	mT(2,1,1)=100.0;
	mT(1,2,1)=100.0;	
	mT(1,1,2)=100.0;
	mT(1,1,1)=100.0;
}

void MACGrid::advectVelocity(double dt)
{
	// TODO: Calculate new velocities and store in target.
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;

	// iterate over each X face
	FOR_EACH_FACE_X {
		// convert cell index to world coordinates
		glm::dvec3 pt = getCenter(i, j, k);
		pt[0] -= theCellSize/2.0f;

		// get full dimensional velocity
		glm::dvec3 vel = getVelocity(pt);

		// do backwards euler step
		glm::dvec3 bes = pt - dt * vel;

		// interpolate new velocity
		glm::dvec3 nvel = getVelocity(bes);

		// store new velocity
		target.mU(i, j, k) = nvel[0];
	}

	// iterate over each Y face
	FOR_EACH_FACE_Y {
		// convert cell index to world coordinates
		glm::dvec3 pt = getCenter(i, j, k);
		pt[1] -= theCellSize/2.0f;

		// get full dimensional velocity
		glm::dvec3 vel = getVelocity(pt);

		// do backwards euler step
		glm::dvec3 bes = pt - dt * vel;

		// interpolate new velocity
		glm::dvec3 nvel = getVelocity(bes);

		// store new velocity
		target.mV(i, j, k) = nvel[1];
	}

	// iterate over each Z face
	FOR_EACH_FACE_Z {
		// convert cell index to world coordinates
		glm::dvec3 pt = getCenter(i, j, k);
		pt[2] -= theCellSize/2.0f;

		// get full dimensional velocity
		glm::dvec3 vel = getVelocity(pt);

		// do backwards euler step
		glm::dvec3 bes = pt - dt * vel;

		// interpolate new velocity
		glm::dvec3 nvel = getVelocity(bes);

		// store new velocity
		target.mW(i, j, k) = nvel[2];
	}

	// Then save the result to our object.
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}

void MACGrid::advectTemperature(double dt)
{
	// TODO: Calculate new temperature and store in target.
	target.mT = mT;

	// temperature is stored per cell
	FOR_EACH_CELL {
		// convert cell index to world coordinates
		glm::dvec3 pt = getCenter(i, j, k);

		// get interpolated velocity
		glm::dvec3 vel = getVelocity(pt);

		// euler step to get previous position
		glm::dvec3 bes = pt - dt * vel;

		// get the interpolated temperature
		target.mT(i,j,k) = getTemperature(bes);
	}

	// Then save the temperature to our object.
	mT = target.mT;
}

void MACGrid::advectDensity(double dt)
{
	// TODO: Calculate new density and store in target.
	target.mD = mD;

	// density is stored per cell
	FOR_EACH_CELL {
		// convert cell index to world coordinates
		glm::dvec3 pt = getCenter(i,j,k);

		// get interpolated velocity at the world point
		glm::dvec3 vel = getVelocity(pt);

		// euler step to get previous position
		glm::dvec3 bes = pt - dt * vel;

		// get the interpolated density
		target.mD(i,j,k) = getDensity(bes);
	}

	// Then save the result to our object.
	mD = target.mD;
}

void MACGrid::computeBuoyancy(double dt)
{
	// TODO: Calculate buoyancy and store in target
	target.mV = mV;

	FOR_EACH_FACE
	{
		if(j==0) continue;

		glm::dvec3 pos = getCenter(i,j,k);
		pos[1] = pos[1] - theCellSize/2.0f;

		double T = getTemperature(pos);
		double D = getDensity(pos);
		double Buoyancy = -theBuoyancyAlpha*D + theBuoyancyBeta*(T - theBuoyancyAmbientTemperature);

		target.mV(i,j,k) = target.getVelocityY(pos) + Buoyancy*dt;
	}
	mV = target.mV;
}

void MACGrid::computeVorticityConfinement(double dt)
{
	// TODO: Calculate vorticity confinement forces.
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;

	GridData omega_i;
	GridData omega_j;
	GridData omega_k;

	omega_i.initialize();
	omega_j.initialize();
	omega_k.initialize();

	FOR_EACH_CELL
	{
		double wjp1 = target.mW(i,j+1,k);
		double wjn1 = target.mW(i,j-1,k);
		double vkp1 = target.mV(i,j,k+1);
		double vkn1 = target.mV(i,j,k-1);
		double ukp1 = target.mU(i,j,k+1);
		double ukn1 = target.mU(i,j,k-1);

		double wip1 = target.mW(i+1,j,k);
		double win1 = target.mW(i-1,j,k);
		double vip1 = target.mV(i+1,j,k);
		double vin1 = target.mV(i-1,j,k);
		double ujp1 = target.mU(i,j+1,k);
		double ujn1 = target.mU(i,j-1,k);

		omega_i(i,j,k) = (wjp1 - wjn1 - vkp1 + vkn1)/(2*theCellSize);
		omega_j(i,j,k) = (ukp1 - ukn1 - wip1 + win1)/(2*theCellSize);
		omega_k(i,j,k) = (vip1 - vin1 - ujp1 + ujn1)/(2*theCellSize);

	}
	GridData F_i;
	GridData F_j;
	GridData F_k;

	F_i.initialize();
	F_j.initialize();
	F_k.initialize();

	FOR_EACH_CELL
	{
		glm::dvec3 wip1(omega_i(i+1,j,k),omega_j(i+1,j,k),omega_k(i+1,j,k));
		glm::dvec3 win1(omega_i(i-1,j,k),omega_j(i-1,j,k),omega_k(i-1,j,k));
		glm::dvec3 wjp1(omega_i(i,j+1,k),omega_j(i,j+1,k),omega_k(i,j+1,k));
		glm::dvec3 wjn1(omega_i(i,j-1,k),omega_j(i,j-1,k),omega_k(i,j-1,k));
		glm::dvec3 wkp1(omega_i(i,j,k+1),omega_j(i,j,k+1),omega_k(i,j,k+1));
		glm::dvec3 wkn1(omega_i(i,j,k-1),omega_j(i,j,k-1),omega_k(i,j,k-1));

		double deta_omega_i = (glm::length(wip1) - glm::length(win1))/(2*theCellSize);
		double deta_omega_j = (glm::length(wjp1) - glm::length(wjn1))/(2*theCellSize);
		double deta_omega_k = (glm::length(wkp1) - glm::length(wkn1))/(2*theCellSize);

		glm::dvec3 deta_omega(deta_omega_i,deta_omega_j,deta_omega_k);
		glm::dvec3 N = deta_omega/(glm::length(deta_omega)+1e-10);
		//std::cout<<N[0]<<" "<<N[1]<<" "<<N[2]<<std::endl;

		glm::dvec3 omega(omega_i(i,j,k),omega_j(i,j,k),omega_k(i,j,k));
		glm::dvec3 F = theVorticityEpsilon*theCellSize*glm::cross(N, omega);

		F_i(i,j,k) = F[0];
		F_j(i,j,k) = F[1];
		F_k(i,j,k) = F[2];
	}

	// Apply the forces to the current velocity and store the result in target.
	FOR_EACH_FACE_X
	{
		if(i==0||i==theDim[MACGrid::X]) {target.mU(i,j,k) = 0;continue;}
		glm::dvec3 p = getCenter(i,j,k);
		p[0] -= theCellSize/2.0f;
		double Fx = F_i.interpolate(p);
		target.mU(i,j,k) = getVelocityX(p) + Fx * dt;
	}

	FOR_EACH_FACE_Y
	{
		if(j==0||j==theDim[MACGrid::Y]) {target.mV(i,j,k) = 0;continue;}
		glm::dvec3 p = getCenter(i,j,k);
		p[1] -= theCellSize/2.0f;
		double Fy = F_i.interpolate(p);
		target.mV(i,j,k) = getVelocityY(p) + Fy * dt;
	}

	FOR_EACH_FACE_Z
	{
		if(k==0||k==theDim[MACGrid::Z]) {target.mW(i,j,k) = 0;continue;}
		glm::dvec3 p = getCenter(i,j,k);
		p[2] -= theCellSize/2.0f;
		double Fz = F_i.interpolate(p);
		target.mW(i,j,k) = getVelocityZ(p) + Fz * dt;
	}

	// Then save the result to our object.
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}

void MACGrid::addExternalForces(double dt)
{
	computeBuoyancy(dt);
	computeVorticityConfinement(dt);
}

void MACGrid::project(double dt)
{
	// TODO: Solve Ax = b for pressure
	// 1. Construct b
	// 2. Construct A 
	// 3. Solve for p
	// Subtract pressure from our velocity and save in target
	// Then save the result to our object

	target.mP = mP;
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;

	// construct b
	GridData b = mP;
	FOR_EACH_CELL {
		// compute constant
		double c = -1.0 * theAirDensity * theCellSize * theCellSize / dt;

		// compute change in x velocity
		double u = (mU(i+1,j,k) - mU(i,j,k)) / theCellSize;
		// compute change in y velocity
		double v = (mV(i,j+1,k) - mV(i,j,k)) / theCellSize;
		// compute change in y velocity
		double w = (mW(i,j,k+1) - mW(i,j,k)) / theCellSize;

		b(i,j,k) = c * (u + v + w);
	}

	// construct A
	//  matrix consisting of pressure neighbor information
	//  already constructed for boundaries
	//  TODO: add support for obstacles in the grid
	GridDataMatrix A = AMatrix;

	// solve for new pressures such that the fluid remains incompressible
	// TODO: what maxIterations and tolerance to use?
	bool div = conjugateGradient(A, target.mP, b, 1e5, 1e-5);
	assert(div==true);

	// update velocities from new pressures
	// vn = v - dt * (1/rho) * dP
	FOR_EACH_CELL {
		// update velocity faces (+1 cell index)
		// boundaries should not change
		// update velocity
		target.mU(i+1,j,k) = mU(i+1,j,k) - dt * (target.mP(i+1,j,k) - target.mP(i,j,k));
		target.mV(i,j+1,k) = mV(i,j+1,k) - dt * (target.mP(i,j+1,k) - target.mP(i,j,k));
		target.mW(i,j,k+1) = mW(i,j,k+1) - dt * (target.mP(i,j,k+1) - target.mP(i,j,k));
	}

	// Then save the result to our object
	mP = target.mP;
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}

glm::dvec3 MACGrid::getVelocity(const glm::dvec3& pt) {
	// TODO : Given a point in space, give the 3D velocity field at the point
	glm::dvec3 vel;
	vel[0] = getVelocityX(pt); 
	vel[1] = getVelocityY(pt); 
	vel[2] = getVelocityZ(pt); 
	return vel;
}

double MACGrid::getVelocityX(const glm::dvec3& pt) {
	return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const glm::dvec3& pt) {
	return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const glm::dvec3& pt) {
	return mW.interpolate(pt);
}

double MACGrid::getTemperature(const glm::dvec3& pt) {
	return mT.interpolate(pt);
}

double MACGrid::getDensity(const glm::dvec3& pt) {
	return mD.interpolate(pt);
}

glm::dvec3 MACGrid::getCenter(int i, int j, int k)
{
	double xstart = theCellSize/2.0;
	double ystart = theCellSize/2.0;
	double zstart = theCellSize/2.0;

	double x = xstart + i*theCellSize;
	double y = ystart + j*theCellSize;
	double z = zstart + k*theCellSize;
	return glm::dvec3(x, y, z);
}

bool MACGrid::isValidCell(int i, int j, int k)
{
	if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
		return false;
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}

void MACGrid::calculateAMatrix() {

	FOR_EACH_CELL {

		int numFluidNeighbors = 0;
		if (i-1 >= 0) {
			AMatrix.plusI(i-1,j,k) = -1;
			numFluidNeighbors++;
		}
		if (i+1 < theDim[MACGrid::X]) {
			AMatrix.plusI(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (j-1 >= 0) {
			AMatrix.plusJ(i,j-1,k) = -1;
			numFluidNeighbors++;
		}
		if (j+1 < theDim[MACGrid::Y]) {
			AMatrix.plusJ(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (k-1 >= 0) {
			AMatrix.plusK(i,j,k-1) = -1;
			numFluidNeighbors++;
		}
		if (k+1 < theDim[MACGrid::Z]) {
			AMatrix.plusK(i,j,k) = -1;
			numFluidNeighbors++;
		}
		// Set the diagonal:
		AMatrix.diag(i,j,k) = numFluidNeighbors;
	}
}

bool MACGrid::conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance) {
	// Solves Ap = d for p.
	FOR_EACH_CELL {
		p(i,j,k) = 0.0; // Initial guess p = 0.	
	}
	GridData r = d; // Residual vector.

	// TODO: Apply a preconditioner here.	
	double tao = 0.97;
	GridData z; z.initialize(0);
	GridData q; q.initialize(0);
	GridData precon; precon.initialize(0);

	FOR_EACH_CELL{
		double ea = A.diag (i,j,k) 
			- (A.plusI (i-1,j,k)*precon (i-1,j,k))*(A.plusI (i-1,j,k)*precon (i-1,j,k))
			- (A.plusJ (i,j-1,k)*precon (i,j-1,k))*(A.plusJ (i,j-1,k)*precon (i,j-1,k))
			- (A.plusK (i,j,k-1)*precon (i,j,k-1))*(A.plusK (i,j,k-1)*precon (i,j,k-1))
			- tao*(A.plusI(i-1,j,k)*(A.plusJ(i-1,j,k) + A.plusK(i-1,j,k))*precon(i-1,j,k)*precon(i-1,j,k)
			+A.plusJ(i,j-1,k)*(A.plusI(i,j-1,k) + A.plusK(i,j-1,k))*precon(i,j-1,k)*precon(i,j-1,k)
			+A.plusK(i,j,k-1)*(A.plusI(i,j,k-1) + A.plusJ(i,j,k-1))*precon(i,j,k-1)*precon(i,j,k-1));
		precon (i,j,k) = 1.0/sqrt(ea+1e-10); 
	}
	//Construct MIC(0) preconditioner
	FOR_EACH_CELL{
		double t = r(i,j,k) 
			- A.plusI(i-1,j,k)*precon(i-1,j,k)*q(i-1,j,k)
			- A.plusJ(i,j-1,k)*precon(i,j-1,k)*q(i,j-1,k)
			- A.plusK(i,j,k-1)*precon(i,j,k-1)*q(i,j,k-1); 
		q(i,j,k) = t*precon(i,j,k); 
	}

	FOR_EACH_CELL_REVERSE{
		double t = q(i,j,k) 
			- A.plusI(i,j,k)*precon(i,j,k)*z(i+1,j,k)
			- A.plusJ(i,j,k)*precon(i,j,k)*z(i,j+1,k)
			- A.plusK(i,j,k)*precon(i,j,k)*z(i,j,k+1); 
		z(i,j,k) = t*precon(i,j,k); 			 
	}

	GridData s = z; // Search vector;
	double sigma = dotProduct(z, r);

	for (int iteration = 0; iteration < maxIterations; iteration++) {

		double rho = sigma;
		apply(A, s, z);//z=applyA(s)
		double alpha = rho/dotProduct(z, s);

		GridData alphaS; alphaS.initialize();
		multiply(alpha, s, alphaS);
		add(p, alphaS, p);
		//p += alpha * s;    

		GridData alphaZ; alphaZ.initialize();
		multiply(alpha, z, alphaZ);
		subtract(r, alphaZ, r);
		//r -= alpha * z;   

		if (maxMagnitude(r) <= tolerance) {
			//PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations."); 
			return true;
		}
		//Construct MIC(0) preconditioner
		FOR_EACH_CELL{
			double t = r(i,j,k) - A.plusI(i-1,j,k)*precon(i-1,j,k)*q(i-1,j,k)
				- A.plusJ(i,j-1,k)*precon(i,j-1,k)*q(i,j-1,k)
				- A.plusK(i,j,k-1)*precon(i,j,k-1)*q(i,j,k-1); 
			q(i,j,k) = t*precon(i,j,k); 
		}

		FOR_EACH_CELL_REVERSE{
			double t = q(i,j,k) 
				- A.plusI(i,j,k)*precon(i,j,k)*z(i+1,j,k)
				- A.plusJ(i,j,k)*precon(i,j,k)*z(i,j+1,k)
				- A.plusK(i,j,k)*precon(i,j,k)*z(i,j,k+1); 
			z(i,j,k) = t*precon(i,j,k); 
		}

		sigma = dotProduct(z, r);
		double beta = sigma / rho;

		GridData betaS; betaS.initialize();
		multiply(beta, s, betaS);
		add(z, betaS, s);
		//s = z + beta * s;    
	}
	PRINT_LINE( "PCG didn't converge!" );
	return false;
}

double MACGrid::dotProduct(const GridData & vector1, const GridData & vector2) {

	double result = 0.0;

	FOR_EACH_CELL {
		result += vector1(i,j,k) * vector2(i,j,k);
	}

	return result;
}

void MACGrid::add(const GridData & vector1, const GridData & vector2, GridData & result) {

	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
	}

}

void MACGrid::subtract(const GridData & vector1, const GridData & vector2, GridData & result) {

	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
	}

}

void MACGrid::multiply(const double scalar, const GridData & vector, GridData & result) {

	FOR_EACH_CELL {
		result(i,j,k) = scalar * vector(i,j,k);
	}

}

double MACGrid::maxMagnitude(const GridData & vector) {

	double result = 0.0;

	FOR_EACH_CELL {
		if (abs(vector(i,j,k)) > result) result = abs(vector(i,j,k));
	}

	return result;
}

void MACGrid::apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result) {

	FOR_EACH_CELL { // For each row of the matrix.

		double diag = 0;
		double plusI = 0;
		double plusJ = 0;
		double plusK = 0;
		double minusI = 0;
		double minusJ = 0;
		double minusK = 0;

		diag = matrix.diag(i,j,k) * vector(i,j,k);
		if (isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
		if (isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
		if (isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
		if (isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
		if (isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
		if (isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);

		result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
	}

}

void MACGrid::saveSmoke(const char* fileName) {
	std::ofstream fileOut(fileName);
	if (fileOut.is_open()) {
		FOR_EACH_CELL {
			fileOut << mD(i,j,k) << std::endl;
		}
		fileOut.close();
	}
}

void MACGrid::draw(const Camera& c)
{   
	drawWireGrid();
	if (theDisplayVel) drawVelocities();   
	if (theRenderMode == CUBES) drawSmokeCubes(c);
	else drawSmoke(c);
}

void MACGrid::drawVelocities()
{
	// draw line at each center
	//glColor4f(0.0, 1.0, 0.0, 1.0);

	glBegin(GL_LINES);
	FOR_EACH_CELL
	{
		glm::dvec3 pos = getCenter(i,j,k);
		glm::dvec3 vel = getVelocity(pos);
		if (glm::length(vel) > 0.0001)
		{
			//vel.Normalize();
			vel *= theCellSize/2.0;
			vel += pos;
			glColor4f(1.0, 1.0, 0.0, 1.0);

			GLdouble doublePos[3];
			doublePos[0] = pos.x, doublePos[1] = pos.y, doublePos[2] = pos.z;
			glVertex3dv(doublePos);

			GLdouble doubleVel[3];
			glColor4f(0.0, 1.0, 0.0, 1.0);
			doubleVel[0] = vel.x, doubleVel[1] = vel.y, doubleVel[2] = vel.z;
			glVertex3dv(doubleVel);
		}
	}
	glEnd();
}

glm::dvec4 MACGrid::getRenderColor(int i, int j, int k)
{
	// TODO : get density value and display as alpha value
	double value = mD(i, j, k); 

#ifdef __TEMP_FIELD__
	double t = mT(i, j, k);
	const double tMax = 100;

	t = std::max(theBuoyancyAmbientTemperature, std::min(t, tMax));
	double frac = (t - theBuoyancyAmbientTemperature)/(tMax - theBuoyancyAmbientTemperature);

	//if(frac > 0.75){	
	//	return glm::dvec4(3 - 3*frac, 4 - 4*frac, frac, value);
	//}else if (frac > 0.5 && frac <= 0.75){
	//	return glm::dvec4(1.5 - frac, 1.0, 2*frac - 0.75, value);
	//}else if(frac > 0.25 && frac <= 0.5){
	//	return glm::dvec4(1.0, 3*frac - 0.5, frac -0.25, value);
	//}else{
	//	return glm::dvec4(4*frac, frac, 0., value);
	//}

	if(frac > 2.0/3){
		return glm::dvec4(1.0, 1.0, 3*frac - 2, value);
	}else if(frac < 1.0/3){
		return glm::dvec4(3*frac, 0.0, 0.0, value);
	}else{
		return glm::dvec4(1.0, 3*frac - 1, 0.0, value);
	}

#else
	return glm::dvec4(1.0, 1.0, 1.0, value);
#endif
}

glm::dvec4 MACGrid::getRenderColor(const glm::dvec3& pt)
{
	// TODO : get density value and display as alpha value
	double value = getDensity(pt); 

#ifdef __TEMP_FIELD__
	double t = getTemperature(pt);
	const double tMax = 100;

	t = std::max(theBuoyancyAmbientTemperature, std::min(t, tMax));
	double frac = (t - theBuoyancyAmbientTemperature)/(tMax - theBuoyancyAmbientTemperature);

	if(frac > 2.0/3){
		return glm::dvec4(1.0, 1.0, 3*frac - 2, value);
	}else if(frac < 1.0/3){
		return glm::dvec4(3*frac, 0.0, 0.0, value);
	}else{
		return glm::dvec4(1.0, 3*frac - 1, 0.0, value);
	}

#else
	return glm::dvec4(1.0, 1.0, 1.0, value);
#endif
}

void MACGrid::drawZSheets(bool backToFront)
{
	// Draw K Sheets from back to front
	double back =  (theDim[2])*theCellSize;
	double top  =  (theDim[1])*theCellSize;
	double right = (theDim[0])*theCellSize;

	double stepsize = theCellSize*0.25;

	double startk = back - stepsize;
	double endk = 0;
	double stepk = -theCellSize;

	if (!backToFront)
	{
		startk = 0;
		endk = back;   
		stepk = theCellSize;
	}

	for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
	{
		for (double j = 0.0; j < top; )
		{
			glBegin(GL_QUAD_STRIP);
			for (double i = 0.0; i <= right; i += stepsize)
			{
				glm::dvec3 pos1 = glm::dvec3(i,j,k); 
				glm::dvec3 pos2 = glm::dvec3(i, j+stepsize, k); 

				glm::dvec4 color1 = getRenderColor(pos1);
				glm::dvec4 color2 = getRenderColor(pos2);

				glColor4dv(glm::value_ptr(color1));
				glVertex3dv(glm::value_ptr(pos1));

				glColor4dv(glm::value_ptr(color2));
				glVertex3dv(glm::value_ptr(pos2));
			} 
			glEnd();
			j+=stepsize;

			glBegin(GL_QUAD_STRIP);
			for (double i = right; i >= 0.0; i -= stepsize)
			{
				glm::dvec3 pos1 = glm::dvec3(i,j,k); 
				glm::dvec3 pos2 = glm::dvec3(i, j+stepsize, k); 

				glm::dvec4 color1 = getRenderColor(pos1);
				glm::dvec4 color2 = getRenderColor(pos2);

				glColor4dv(glm::value_ptr(color1));
				glVertex3dv(glm::value_ptr(pos1));

				glColor4dv(glm::value_ptr(color2));
				glVertex3dv(glm::value_ptr(pos2));
			} 
			glEnd();
			j+=stepsize;
		}
	}
}

void MACGrid::drawXSheets(bool backToFront)
{
	// Draw K Sheets from back to front
	double back =  (theDim[2])*theCellSize;
	double top  =  (theDim[1])*theCellSize;
	double right = (theDim[0])*theCellSize;

	double stepsize = theCellSize*0.25;

	double starti = right - stepsize;
	double endi = 0;
	double stepi = -theCellSize;

	if (!backToFront)
	{
		starti = 0;
		endi = right;   
		stepi = theCellSize;
	}

	for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
	{
		for (double j = 0.0; j < top; )
		{
			glBegin(GL_QUAD_STRIP);
			for (double k = 0.0; k <= back; k += stepsize)
			{
				glm::dvec3 pos1 = glm::dvec3(i,j,k); 
				glm::dvec3 pos2 = glm::dvec3(i, j+stepsize, k); 

				glm::dvec4 color1 = getRenderColor(pos1);
				glm::dvec4 color2 = getRenderColor(pos2);

				glColor4dv(glm::value_ptr(color1));
				glVertex3dv(glm::value_ptr(pos1));

				glColor4dv(glm::value_ptr(color2));
				glVertex3dv(glm::value_ptr(pos2));
			} 
			glEnd();
			j+=stepsize;

			glBegin(GL_QUAD_STRIP);
			for (double k = back; k >= 0.0; k -= stepsize)
			{
				glm::dvec3 pos1 = glm::dvec3(i,j,k); 
				glm::dvec3 pos2 = glm::dvec3(i, j+stepsize, k); 

				glm::dvec4 color1 = getRenderColor(pos1);
				glm::dvec4 color2 = getRenderColor(pos2);

				glColor4dv(glm::value_ptr(color1));
				glVertex3dv(glm::value_ptr(pos1));

				glColor4dv(glm::value_ptr(color2));
				glVertex3dv(glm::value_ptr(pos2));
			} 
			glEnd();
			j+=stepsize;
		}
	}
}


void MACGrid::drawSmoke(const Camera& c)
{
	glm::dvec3 eyeDir = c.getBackward();
	double zresult = fabs(glm::dot(eyeDir, glm::dvec3(1,0,0)));
	double xresult = fabs(glm::dot(eyeDir, glm::dvec3(0,0,1)));
	//double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

	if (zresult < xresult)
	{      
		drawZSheets(c.getPosition()[2] < 0);
	}
	else 
	{
		drawXSheets(c.getPosition()[0] < 0);
	}
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
	std::multimap<double, MACGrid::Cube, std::greater<double> > cubes;
	FOR_EACH_CELL
	{
		MACGrid::Cube cube;
		cube.color = getRenderColor(i,j,k);
		cube.pos = getCenter(i,j,k);
		cube.dist = glm::length((cube.pos - c.getPosition()));
		cubes.insert(std::make_pair(cube.dist, cube));
	} 

	// Draw cubes from back to front
	std::multimap<double, MACGrid::Cube, std::greater<double> >::const_iterator it;
	for (it = cubes.begin(); it != cubes.end(); ++it)
	{
		drawCube(it->second);
	}
}

void MACGrid::drawWireGrid()
{
	// Display grid in light grey, draw top & bottom

	double xstart = 0.0;
	double ystart = 0.0;
	double zstart = 0.0;
	double xend = theDim[0]*theCellSize;
	double yend = theDim[1]*theCellSize;
	double zend = theDim[2]*theCellSize;

	glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
	glDisable(GL_LIGHTING);
	glColor3f(0.25, 0.25, 0.25);

	glBegin(GL_LINES);
	for (int i = 0; i <= theDim[0]; i++)
	{
		double x = xstart + i*theCellSize;
		glVertex3d(x, ystart, zstart);
		glVertex3d(x, ystart, zend);

		glVertex3d(x, yend, zstart);
		glVertex3d(x, yend, zend);
	}

	for (int i = 0; i <= theDim[2]; i++)
	{
		double z = zstart + i*theCellSize;
		glVertex3d(xstart, ystart, z);
		glVertex3d(xend, ystart, z);

		glVertex3d(xstart, yend, z);
		glVertex3d(xend, yend, z);
	}

	glVertex3d(xstart, ystart, zstart);
	glVertex3d(xstart, yend, zstart);

	glVertex3d(xend, ystart, zstart);
	glVertex3d(xend, yend, zstart);

	glVertex3d(xstart, ystart, zend);
	glVertex3d(xstart, yend, zend);

	glVertex3d(xend, ystart, zend);
	glVertex3d(xend, yend, zend);
	glEnd();
	glPopAttrib();

	glEnd();
}

#define LEN 0.5
void MACGrid::drawFace(const MACGrid::Cube& cube)
{
	glColor4dv(glm::value_ptr(cube.color));
	glPushMatrix();
	glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
	glScaled(theCellSize, theCellSize, theCellSize);
	glBegin(GL_QUADS);
	glNormal3d( 0.0,  0.0, 1.0);
	glVertex3d(-LEN, -LEN, LEN);
	glVertex3d(-LEN,  LEN, LEN);
	glVertex3d( LEN,  LEN, LEN);
	glVertex3d( LEN, -LEN, LEN);
	glEnd();
	glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube& cube)
{
	glColor4dv(glm::value_ptr(cube.color));
	glPushMatrix();
	glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
	glScaled(theCellSize, theCellSize, theCellSize);
	glBegin(GL_QUADS);
	glNormal3d( 0.0, -1.0,  0.0);
	glVertex3d(-LEN, -LEN, -LEN);
	glVertex3d(-LEN, -LEN,  LEN);
	glVertex3d( LEN, -LEN,  LEN);
	glVertex3d( LEN, -LEN, -LEN);         

	glNormal3d( 0.0,  0.0, -0.0);
	glVertex3d(-LEN, -LEN, -LEN);
	glVertex3d(-LEN,  LEN, -LEN);
	glVertex3d( LEN,  LEN, -LEN);
	glVertex3d( LEN, -LEN, -LEN);

	glNormal3d(-1.0,  0.0,  0.0);
	glVertex3d(-LEN, -LEN, -LEN);
	glVertex3d(-LEN, -LEN,  LEN);
	glVertex3d(-LEN,  LEN,  LEN);
	glVertex3d(-LEN,  LEN, -LEN);

	glNormal3d( 0.0, 1.0,  0.0);
	glVertex3d(-LEN, LEN, -LEN);
	glVertex3d(-LEN, LEN,  LEN);
	glVertex3d( LEN, LEN,  LEN);
	glVertex3d( LEN, LEN, -LEN);

	glNormal3d( 0.0,  0.0, 1.0);
	glVertex3d(-LEN, -LEN, LEN);
	glVertex3d(-LEN,  LEN, LEN);
	glVertex3d( LEN,  LEN, LEN);
	glVertex3d( LEN, -LEN, LEN);

	glNormal3d(1.0,  0.0,  0.0);
	glVertex3d(LEN, -LEN, -LEN);
	glVertex3d(LEN, -LEN,  LEN);
	glVertex3d(LEN,  LEN,  LEN);
	glVertex3d(LEN,  LEN, -LEN);
	glEnd();
	glPopMatrix();
}