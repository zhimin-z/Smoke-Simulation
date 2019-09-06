#include "grid_data.h"
//#define __CUBIC_INTERP__

GridData::GridData() :
	mDfltValue(0.0), mMax(0.0,0.0,0.0){
		// TODO : GridData constructor
}

GridData::GridData(const GridData& orig) {
	// TODO : GridData copy constructor
	mData = orig.mData;
	mMax = orig.mMax;
}

GridData::~GridData() {
	// TODO : GridData destructor
}

std::vector<double>& GridData::data() {
	// TODO : Return underlying data structure (you may change the method header
	// to fit whatever design you choose).
	return mData;
}

GridData& GridData::operator=(const GridData& orig) {
	// TODO : Override GridData '=' operator with copy 
	if (this != &orig){
		mDfltValue = orig.mDfltValue;
		mData = orig.mData;
		mMax = orig.mMax;
	}
	return *this;
}

void GridData::initialize(double dfltValue) {
	// TODO : Initialize the grid to a default value
	mDfltValue = dfltValue;
	mMax[0] = theCellSize*theDim[0];
	mMax[1] = theCellSize*theDim[1];
	mMax[2] = theCellSize*theDim[2];
	mData.resize(theDim[0]*theDim[1]*theDim[2], false);
	std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridData::operator()(int i, int j, int k) {
	// TODO : Grid accessor that allows for client to access and set cell data
	static double dflt = mDfltValue; 
	if (i< 0 || j<0 || k<0 || 
		i >= theDim[0] || 
		j >= theDim[1] || 
		k >= theDim[2]) 
		return dflt;

	int col = i;
	int row = k*theDim[0];
	int stack = j*theDim[0]*theDim[2];

	return mData[col+row+stack];
}

const double GridData::operator()(int i, int j, int k) const {
	// TODO : Grid accessor
	static double dflt = mDfltValue; 
	if (i< 0 || j<0 || k<0 || 
		i >= theDim[0] || 
		j >= theDim[1] || 
		k >= theDim[2]) 
		return dflt;

	int col = i;
	int row = k*theDim[0];
	int stack = j*theDim[0]*theDim[2];

	return mData[col+row+stack];
}

void GridData::getCell(const glm::dvec3& pt, int& i, int& j, int& k) {
	// TODO : Given a point in world coordinates, return the cell index
	// corresponding to it.
	i = (int) (pt[0]/theCellSize);
	j = (int) (pt[1]/theCellSize);
	k = (int) (pt[2]/theCellSize); 
}

double GridData::interpolate(const glm::dvec3& pt) {
	// TODO : Given a point, interpolate the value in the grid at that point
	int i, j, k;
	glm::dvec3 pos = worldToSelf(pt);
	getCell(pos, i, j, k);

	glm::dvec3 frac;
	double scale = 1.0/theCellSize;  
	frac[0] = scale*(pos[0] - i*theCellSize);
	frac[1] = scale*(pos[1] - j*theCellSize);
	frac[2] = scale*(pos[2] - k*theCellSize);

	//std::cout<<frac[0]<<" "<<frac[1]<<" "<<frac[2]<<std::endl;

	assert (frac[0] < 1.0 && frac[0] >= 0);
	assert (frac[1] < 1.0 && frac[1] >= 0);	
	assert (frac[2] < 1.0 && frac[2] >= 0);
	
#ifdef __CUBIC_INTERP__
	// Y interpolation
	// near-right column
	double fm1m1m1  = (*this)(i-1,  j-1,  k-1);
	double fm10m1   = (*this)(i-1,  j,  k-1);
	double fm11m1   = (*this)(i-1,  j+1,  k-1);
	double fm12m1   = (*this)(i-1,  j+2,  k-1);
	double yi_m10m1  = HCInterpolate(fm1m1m1, fm10m1, fm11m1, fm12m1, frac[1]);
	// near-center column
	double f0m1m1 = (*this)(i,  j-1,  k-1);
	double f00m1  = (*this)(i,    j,  k-1);
	double f01m1  = (*this)(i,  j+1,  k-1);
	double f02m1  = (*this)(i,  j+2,  k-1);
	double yi_00m1 = HCInterpolate(f0m1m1, f00m1, f01m1, f02m1, frac[1]);
	// near-left column
	double f1m1m1 = (*this)(i+1,  j-1,  k-1);
	double f10m1  = (*this)(i+1,    j,  k-1);
	double f11m1  = (*this)(i+1,  j+1,  k-1);
	double f12m1  = (*this)(i+1,  j+2,  k-1);
	double yi_10m1 = HCInterpolate(f1m1m1, f10m1, f11m1, f12m1, frac[1]);
	// near-2left column
	double f2m1m1 = (*this)(i+2,  j-1,  k-1);
	double f20m1  = (*this)(i+2,    j,  k-1);
	double f21m1  = (*this)(i+2,  j+1,  k-1);
	double f22m1  = (*this)(i+2,  j+2,  k-1);
	double yi_20m1 = HCInterpolate(f2m1m1, f20m1, f21m1, f22m1, frac[1]);

	// center-right column
	double fm1m10 = (*this)(i-1,  j-1,  k);
	double fm100  = (*this)(i-1,    j,  k);
	double fm110  = (*this)(i-1,  j+1,  k);
	double fm120  = (*this)(i-1,  j+2,  k);
	double yi_m100 = HCInterpolate(fm1m10, fm100, fm110, fm120, frac[1]);
	// center column
	double f0m10  = (*this)(i,  j-1,  k);
	double f000   = (*this)(i,    j,  k);
	double f010   = (*this)(i,  j+1,  k);
	double f020   = (*this)(i,  j+2,  k);
	double yi_000  = HCInterpolate(f0m10, f000, f010, f020, frac[1]);
	// center-left column
	double f1m10  = (*this)(i+1,  j-1,  k);
	double f100   = (*this)(i+1,    j,  k);
	double f110   = (*this)(i+1,  j+1,  k);
	double f120   = (*this)(i+1,  j+2,  k);
	double yi_100  = HCInterpolate(f1m10, f100, f110, f120, frac[1]);
	// center-2left column
	double f2m10  = (*this)(i+2,  j-1,  k);
	double f200   = (*this)(i+2,    j,  k);
	double f210   = (*this)(i+2,  j+1,  k);
	double f220   = (*this)(i+2,  j+2,  k);
	double yi_200  = HCInterpolate(f2m10, f200, f210, f220, frac[1]);

	// far-right column
	double fm1m11 = (*this)(i-1,  j-1,  k+1);
	double fm101  = (*this)(i-1,    j,  k+1);
	double fm111  = (*this)(i-1,  j+1,  k+1);
	double fm121  = (*this)(i-1,  j+2,  k+1);
	double yi_m101 = HCInterpolate(fm1m11, fm101, fm111, fm121, frac[1]);
	// far-center column
	double f0m11  = (*this)(i,  j-1,  k+1);
	double f001   = (*this)(i,    j,  k+1);
	double f011   = (*this)(i,  j+1,  k+1);
	double f021   = (*this)(i,  j+2,  k+1);
	double yi_001  = HCInterpolate(f0m11, f001, f011, f021, frac[1]);
	// far-left column
	double f1m11  = (*this)(i+1,  j-1,  k+1);
	double f101   = (*this)(i+1,    j,  k+1);
	double f111   = (*this)(i+1,  j+1,  k+1);
	double f121   = (*this)(i+1,  j+2,  k+1);
	double yi_101  = HCInterpolate(f1m11, f101, f111, f121, frac[1]);
	// far-2left column
	double f2m11  = (*this)(i+2,  j-1,  k+1);
	double f201   = (*this)(i+2,    j,  k+1);
	double f211   = (*this)(i+2,  j+1,  k+1);
	double f221   = (*this)(i+2,  j+2,  k+1);
	double yi_201  = HCInterpolate(f2m11, f201, f211, f221, frac[1]);

	// +2 far-right column
	double fm1m12   = (*this)(i-1,  j-1,  k+2);
	double fm102    = (*this)(i-1,    j,  k+2);
	double fm112    = (*this)(i-1,  j+1,  k+2);
	double fm122    = (*this)(i-1,  j+2,  k+2);
	double yi_m102   = HCInterpolate(fm1m12, fm102, fm112, fm122, frac[1]);
	// +2 far-center column
	double f0m12    = (*this)(i,  j-1,  k+2);
	double f002     = (*this)(i,    j,  k+2);
	double f012     = (*this)(i,  j+1,  k+2);
	double f022     = (*this)(i,  j+2,  k+2);
	double yi_002    = HCInterpolate(f0m12, f002, f012, f022, frac[1]);
	// +2 far-left column
	double f1m12    = (*this)(i+1,  j-1,  k+2);
	double f102     = (*this)(i+1,    j,  k+2);
	double f112     = (*this)(i+1,  j+1,  k+2);
	double f122     = (*this)(i+1,  j+2,  k+2);
	double yi_102    = HCInterpolate(f1m12, f102, f112, f122, frac[1]);
	// +2 far-2left column
	double f2m12    = (*this)(i+2,  j-1,  k+2);
	double f202     = (*this)(i+2,    j,  k+2);
	double f212     = (*this)(i+2,  j+1,  k+2);
	double f222     = (*this)(i+2,  j+2,  k+2);
	double yi_202    = HCInterpolate(f2m12, f202, f212, f222, frac[1]);

	// Z interpolation
	// right column
	double zi_m100 = HCInterpolate(yi_m10m1, yi_m100, yi_m101, yi_m102, frac[2]);
	// center column
	double zi_000 = HCInterpolate(yi_00m1, yi_000, yi_001, yi_002, frac[2]);
	// left column
	double zi_100 = HCInterpolate(yi_10m1, yi_100, yi_101, yi_102, frac[2]);
	// +2 left column
	double zi_200 = HCInterpolate(yi_20m1, yi_200, yi_201, yi_202, frac[2]);

	// X interpolation
	double xi = HCInterpolate(zi_m100, zi_000, zi_100, zi_200, frac[0]);
	return xi;

	  #else
	double tmp000 = (*this)(i, j, k);
	double tmp100 = (*this)(i+1, j, k);
	double tmp010 = (*this)(i, j+1, k);
	double tmp110 = (*this)(i+1, j+1, k);
	double tmp001 = (*this)(i, j, k+1);
	double tmp011 = (*this)(i, j+1, k+1);
	double tmp101 = (*this)(i+1, j, k+1);
	double tmp111 = (*this)(i+1, j+1, k+1);

	//Interpolate X
	double tmpX00 = LERP(tmp000, tmp100, frac[0]);
	double tmpX10 = LERP(tmp010, tmp110, frac[0]);
	double tmpX01 = LERP(tmp001, tmp101, frac[0]);
	double tmpX11 = LERP(tmp011, tmp111, frac[0]);

	//Interpolate Y
	double tmpXY0 = LERP (tmpX00, tmpX10, frac[1]);
	double tmpXY1 = LERP (tmpX01, tmpX11, frac[1]);
		
	//Interpolate Z
	double tmp= LERP(tmpXY0, tmpXY1, frac[2]);
	return tmp;
	  #endif
}

glm::dvec3 GridData::worldToSelf(const glm::dvec3& pt) const {
	// TODO : Given a point, returns the cell index that the grid uses in its own
	// space.
	glm::dvec3 pos;
	pos[0] = min(max(0.0, pt[0] - theCellSize*0.5), mMax[0]);
	pos[1] = min(max(0.0, pt[1] - theCellSize*0.5), mMax[1]);
	pos[2] = min(max(0.0, pt[2] - theCellSize*0.5), mMax[2]);
	return pos;
}

GridDataX::GridDataX() : GridData() {
}

GridDataX::~GridDataX() {
}

void GridDataX::initialize(double dfltValue) {
	// TODO : Initialize GridDataX
	mDfltValue = dfltValue;
	mMax[0] = theCellSize*(theDim[0]+1);
	mMax[1] = theCellSize*theDim[1];
	mMax[2] = theCellSize*theDim[2];
	mData.resize((theDim[0]+1)*theDim[1]*theDim[2], false);
	std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataX::operator()(int i, int j, int k) {
	// TODO : GridX accessor
	static double dflt = mDfltValue; 

	if (i< 0 || i > theDim[0]) return dflt;

	if (j < 0) j = 0;
	else if (j >= theDim[1]) j = theDim[1]-1;

	if (k < 0) k = 0;
	else if (k >= theDim[2]) k = theDim[2]-1;

	int col = i;
	int row = k*(theDim[0]+1);
	int stack = j*(theDim[0]+1)*theDim[2];

	return mData[col+row+stack];
}

const double GridDataX::operator()(int i, int j, int k) const {
	// TODO : GridX accessor
	static double dflt = mDfltValue; 

	if (i< 0 || i > theDim[0]) return dflt;

	if (j < 0) j = 0;
	else if (j >= theDim[1]) j = theDim[1]-1;

	if (k < 0) k = 0;
	else if (k >= theDim[2]) k = theDim[2]-1;

	int col = i;
	int row = k*(theDim[0]+1);
	int stack = j*(theDim[0]+1)*theDim[2];

	return mData[col+row+stack];
}

glm::dvec3 GridDataX::worldToSelf(const glm::dvec3& pt) const {
	// TODO : Given a point, returns the cell index that the grid uses in its own
	// space
	glm::dvec3 pos;
	pos[0] = min(max(0.0, pt[0]), mMax[0]);
	pos[1] = min(max(0.0, pt[1] - theCellSize*0.5), mMax[1]);
	pos[2] = min(max(0.0, pt[2] - theCellSize*0.5), mMax[2]);
	return pos;
}

GridDataY::GridDataY() : GridData() {
}

GridDataY::~GridDataY() {
}

void GridDataY::initialize(double dfltValue) {
	// TODO : Initialize GridDataY
	mDfltValue = dfltValue;
	mMax[0] = theCellSize*theDim[0];
	mMax[1] = theCellSize*(theDim[1]+1);
	mMax[2] = theCellSize*theDim[2];
	mData.resize(theDim[0]*(theDim[1]+1)*theDim[2], false);
	std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataY::operator()(int i, int j, int k) {
	// TODO : GridY accessor
	static double dflt = mDfltValue; 

	if (j< 0 || j > theDim[1]) return dflt;

	if (i < 0) i = 0;
	else if (i >= theDim[0]) i = theDim[0]-1;

	if (k < 0) k = 0;
	else if (k >= theDim[2]) k = theDim[2]-1;

	int col = i;
	int row = k*theDim[0];
	int stack = j*theDim[0]*theDim[2];

	return mData[col+row+stack];
}

const double GridDataY::operator()(int i, int j, int k) const {
	// TODO : GridY accessor
	static double dflt = mDfltValue; 

	if (j< 0 || j > theDim[1]) return dflt;

	if (i < 0) i = 0;
	else if (i >= theDim[0]) i = theDim[0]-1;

	if (k < 0) k = 0;
	else if (k >= theDim[2]) k = theDim[2]-1;

	int col = i;
	int row = k*theDim[0];
	int stack = j*theDim[0]*theDim[2];

	return mData[col+row+stack];
}

glm::dvec3 GridDataY::worldToSelf(const glm::dvec3& pt) const {
	// TODO : Given a point, returns the cell index that the grid uses in its own
	// space
	glm::dvec3 pos;
	pos[0] = min(max(0.0, pt[0] - theCellSize*0.5), mMax[0]);
	pos[1] = min(max(0.0, pt[1]), mMax[1]);
	pos[2] = min(max(0.0, pt[2] - theCellSize*0.5), mMax[2]);
	return pos;
}

GridDataZ::GridDataZ() : GridData() {
}

GridDataZ::~GridDataZ() {
}

void GridDataZ::initialize(double dfltValue) {
	// TODO : Intialize GridDataZ
	mDfltValue = dfltValue;
	mMax[0] = theCellSize*theDim[0];
	mMax[1] = theCellSize*theDim[1];
	mMax[2] = theCellSize*(theDim[2]+1);
	mData.resize(theDim[0]*theDim[1]*(theDim[2]+1), false);
	std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataZ::operator()(int i, int j, int k) {
	// TODO : GridZ accessor
	static double dflt = mDfltValue; 

	if (k< 0 || k > theDim[2]) return dflt;

	if (i < 0) i = 0;
	else if (i >= theDim[0]) i = theDim[0]-1;

	if (j < 0) j = 0;
	else if (j >= theDim[1]) j = theDim[1]-1;

	int col = i;
	int row = k*theDim[0];
	int stack = j*theDim[0]*(theDim[2]+1);

	return mData[col+row+stack];
}

const double GridDataZ::operator()(int i, int j, int k) const {
	// TODO : GridY accessor
	static double dflt = mDfltValue; 

	if (k< 0 || k > theDim[2]) return dflt;

	if (i < 0) i = 0;
	else if (i >= theDim[0]) i = theDim[0]-1;

	if (j < 0) j = 0;
	else if (j >= theDim[1]) j = theDim[1]-1;

	int col = i;
	int row = k*theDim[0];
	int stack = j*theDim[0]*(theDim[2]+1);

	return mData[col+row+stack];
}

glm::dvec3 GridDataZ::worldToSelf(const glm::dvec3& pt) const {
	// TODO : Given a point, returns the cell index that the grid uses in its own
	// space
	glm::dvec3 pos;
	pos[0] = min(max(0.0, pt[0] - theCellSize*0.5), mMax[0]);
	pos[1] = min(max(0.0, pt[1] - theCellSize*0.5), mMax[1]);
	pos[2] = min(max(0.0, pt[2]), mMax[2]);
	return pos;
}

double GridData::HCInterpolate(double q_1, double q0, double q1, double q2, double x) {
	//Hermite cubic interpolation
	double deltaq = q1 - q0;
	double d0 = (q1 - q_1)/2.0;
	double d1 = (q2 - q0)/2.0;

	if (SIGN(d0) != SIGN(deltaq)) {
		d0 = 0.0;
	}
	if (SIGN(d1) != SIGN(deltaq)) {
		d1 = 0.0;
	}

	double a0 = q0;
	double a1 = d0;
	double a2 = 3*deltaq - 2*d0 - d1;
	double a3 = -2*deltaq + d0 + d1;

	double tmp = a3*x*x*x + a2*x*x + a1*x + a0;
	return tmp;
}