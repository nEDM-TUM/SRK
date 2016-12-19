////////////////////////////////////////////////////////////////
/// class SRKROOTField
/// Contains a class that defines a vector/scalar field in 2D cylindrical symmetry or in 3D
/// and can either linearly or cubically interpolate
/// Requires that the ROOT library be installed and linked
///
/// Author: Matthew Bales
///////////////////////////////////////////////////////////////

#ifndef SRKROOTField_HH
#define SRKROOTField_HH

//Standard Libraries
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/stat.h>

//Root libraries
#include "TROOT.h"
#include "TH3.h"
#include "TH2.h"
#include "TH1.h"
#include "TVector3.h"
#include "TFile.h"
#include "TSystem.h"

#include "SRKODEState.h"

class SRKROOTField
{

public:
	SRKROOTField();
	SRKROOTField(std::string filePath, std::string histName, double scalingValue, int inpSpaceDim, int inpFieldDim);
	~SRKROOTField();
	void reset();
	void linearInterp3D(TVector3 pos, TVector3& vecOut) const;
	void cubicInterp3D(TVector3 pos, TVector3& vecOut) const;
	int loadFieldFromFile(std::string filePath, std::string histName, double scalingValue, int inpSpaceDim, int inpFieldDim);
	bool isPositionInsideField(TVector3 pos);
	int saveFieldToFile(std::string filePath, std::string histName);
	bool fieldFileExists(std::string strFileName);

	//From http://www.paulinternet.nl/?page=bicubic
	TVector3 cubicInterpolate(TVector3 p[4], double xIn) const;  //cubic interpolation along 1D, p is an array of data values around the position, xIn is the value from 0,1 representing how far betweeen p[1] and p[2] we are interpolating
	TVector3 bicubicInterpolate(TVector3 p[4][4], double xIn, double yIn) const; //cubic interpolation in 2D
	TVector3 tricubicInterpolate(TVector3 p[4][4][4], double xIn, double yIn, double zIn) const; //cubic interpolation in 3D


    inline int getRows(){return numBins[0];};
    inline int getColumns(){return numBins[1];};
    inline int getLayers(){return numBins[2];};
    inline double getXStart(){return startCorner.X();};
    inline double getYStart(){return startCorner.Y();};
    inline double getZStart(){return startCorner.Z();};
    inline double getXSpacing(){return spacing.X();};
    inline double getYSpacing(){return spacing.Y();};
    inline double getZSpacing(){return spacing.Z();};
    inline double getXFinal(){return endCorner.X();};
    inline double getYFinal(){return endCorner.Y();};
    inline double getZFinal(){return endCorner.Z();};
    inline double getXLength(){return length.X();};
    inline double getYLength(){return length.Y();};
    inline double getZLength(){return length.Z();};
    inline void setSymmetry(bool inpX, bool inpY=false, bool inpZ=false){symmetry[0]=inpX;symmetry[1]=inpY;symmetry[2]=inpZ; }
    inline void setSymmetry(int symmetryToSet, bool inpBool){symmetry[symmetryToSet]=inpBool;};
    inline void setRotation(double inpx,double inpy,double inpz){rotateX=inpx; rotateY=inpy; rotateZ=inpz;};
    inline TVector3 getVector(int x, int y, int z) const {return TVector3(theField[0][x][y][z],theField[1][x][y][z],theField[2][x][y][z]);};

private:
	int numBins[3];
	TVector3 startCorner;
	TVector3 endCorner;
	TVector3 spacing;
	TVector3 length;

	bool fieldLoaded;
	bool symmetry[3];  //Are we mirroring across any of the axis.  Presumes start corner is at 0 in that axis and that it is mirrored in the negative direction

	std::vector<std::vector<std::vector<double> > > theField[3]; //The field.  Stored as nested STL vectors.  Access by theField[i][x][y][z].

	int dimOfSpace; //2 = cylindrical coordinates(r,z)   3 = 3D cartesian (x,y,z)
	int dimOfField; //Scalar field or vector field?

	//For interpolation
//	double radius;
//    TVector3 resultVec;
//    double modx,mody,modz,modxm,modym,modzm;
//    double dx,dy,dz;
//    int x,y,z,xp1,yp1,zp1;
//    double a1,a2,a3,a4,a5,a6,a7,a8;
	double rotateX, rotateY, rotateZ;  //Rotate angle theta aroudn each axis in order to transform coordinates.  Transformed in order X,Y,Z.

	//private functions
	int loadFieldROOTFile(std::string filePath, std::string histName, double scalingValue);
	int loadFieldTXTFile(std::string filePath, double scalingValue);
	TVector3 getPosInField(TVector3 inp);  //Rotates the mirrors position as well as turns 3D position into radius and z coordinate (only first 2 coordinates used) in the case of 2D space field
	void setSize(int, int, int);

};



#endif // SRKROOTField_HH
