#ifndef MACGrid_H_
#define MACGrid_H_

#pragma warning(disable: 4244 4267 4996)

#include "open_gl_headers.h" // PETER KUTZ.

#include <glm/glm.hpp>
#include <functional>

#include "grid_data.h"
#include "grid_data_matrix.h" // PETER KUTZ.

class Camera;

class MACGrid
{

public:
	MACGrid();
	~MACGrid();
	MACGrid(const MACGrid& orig);
	MACGrid& operator=(const MACGrid& orig);

	void reset();

	void draw(const Camera& c);
	void updateSources();
	void advectVelocity(double dt);
	void addExternalForces(double dt);
	void project(double dt);
	void advectTemperature(double dt);
	void advectDensity(double dt);

protected:

	// Setup
	void initialize();

	// Simulation
	void computeBuoyancy(double dt);
	void computeVorticityConfinement(double dt);

	// Rendering
	struct Cube { glm::dvec3 pos; glm::dvec4 color; double dist; };
	void drawWireGrid();
	void drawSmokeCubes(const Camera& c);
	void drawSmoke(const Camera& c);
	void drawCube(const MACGrid::Cube& c);
	void drawFace(const MACGrid::Cube& c);
	void drawVelocities();
	glm::dvec4 getRenderColor(int i, int j, int k);
	glm::dvec4 getRenderColor(const glm::dvec3& pt);
	void drawZSheets(bool backToFront);
	void drawXSheets(bool backToFront);

	// GridData accessors
	enum Direction { X, Y, Z };
	double getTemperature(const glm::dvec3& pt);
	double getDensity(const glm::dvec3& pt);
	double getVelocityX(const glm::dvec3& pt);
	double getVelocityY(const glm::dvec3& pt);
	double getVelocityZ(const glm::dvec3& pt);
    glm::dvec3 getVelocity(const glm::dvec3& pt);
	glm::dvec3 getCenter(int i, int j, int k);
    bool isValidCell(int i, int j, int k);

  // Sets up A matrix for calculation
	void calculateAMatrix();

  // Conjugate Gradient stuff
	bool conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance);
	double dotProduct(const GridData & vector1, const GridData & vector2);
	void add(const GridData & vector1, const GridData & vector2, GridData & result);
	void subtract(const GridData & vector1, const GridData & vector2, GridData & result);
	void multiply(const double scalar, const GridData & vector, GridData & result);
	double maxMagnitude(const GridData & vector);
	void apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result);

  // TODO : Fill in the necessary data structures to maintain velocity, pressure
  // and density
	GridDataX mU; // X component of velocity, stored on X faces, size is (dimX+1)*dimY*dimZ
	GridDataY mV; // Y component of velocity, stored on Y faces, size is dimX*(dimY+1)*dimZ
	GridDataZ mW; // W component of velocity, stored on Z faces, size is dimX*dimY*(dimZ+1)
	GridData mP;  // Pressure, stored at grid centers, size is dimX*dimY*dimZ
	GridData mD;  // Density, stored at grid centers, size is dimX*dimY*dimZ
	GridData mT;  // Temperature, stored at grid centers, size is dimX*dimY*dimZ

	GridDataMatrix AMatrix;

public:

	enum RenderMode { CUBES, SHEETS };
	static RenderMode theRenderMode;
	static bool theDisplayVel;
	
	void saveSmoke(const char* fileName);
};

#endif
