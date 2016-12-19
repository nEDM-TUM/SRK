#ifndef SRKField_h
#define SRKField_h 1

#include "SRKROOTField.h"

#include <string>

const static double inv_c_light = 1. / 2.99792458e+8;

enum FieldType{FIELD_MAGNETIC,FIELD_ELECTRIC,FIELD_GRAVITY};
enum FieldClass{FIELDCLASS_UNIFORM,FIELDCLASS_GRADIENT,FIELDCLASS_DIPOLE,FIELDCLASS_OSCILLATING,FIELDCLASS_INTERPOLATION,FIELDCLASS_QUADRUPOLE,FIELDCLASS_SEXTUPOLE};

////////////////////////////////////////////////////////////////
/// SRKFieldSettings
/// Simple class to store settings related to SRKFields
///
/// Author: Matthew Bales
///////////////////////////////////////////////////////////////
class SRKFieldSettings
{
public:

	inline TString getFieldTypeString()
	{
		TString out = "Magnetic";
		if(fieldType == FIELD_ELECTRIC) out = "Electric";
		if(fieldType == FIELD_GRAVITY) out = "Gravity";
		return out;
	}
	inline TString getFieldClassString()
	{
		TString out = "Uniform";
		if(fieldClass == FIELDCLASS_GRADIENT) out = "Gradient";
		if(fieldClass == FIELDCLASS_DIPOLE) out = "Dipole";
		if(fieldClass == FIELDCLASS_OSCILLATING) out = "Oscillating";
		if(fieldClass == FIELDCLASS_INTERPOLATION) out = "Interpolation";
		if(fieldClass == FIELDCLASS_QUADRUPOLE) out = "Quadrupole";
		if(fieldClass == FIELDCLASS_SEXTUPOLE) out = "Sextupole";
		return out;
	}
	std::string fieldFilePath; //Path for field file (e.g. for interpolated)
	std::string histName = "field"; //In case of stored histogram, the histogram name

	TVector3 centerPos;  //Center position of field extents. Base for local unit system local unit system for field. Not used by interpolated field.[mm]
	TVector3 extents = TVector3(100., 100., 100.);    //Extents of the field in both directions from the center position (later to be changed)
	TVector3 offset;     //Offset of local unit system
	TVector3 direction;  //Direction of field (changed to unit vector)
	TVector3 axisDirection; //Direction of axis for radial-symmetric fields, e.g. sextupole
	TVector3 moment;	  //

	double scalingValue = 0.; //[tesla][V/m][m/s2] Zero scaling value will be the signal that no field should be loaded;
	double frequency = 0.; 		  //rad/s
	double angleX=0;
	double angleY=0;
	double angleZ=0;
	int spaceDim = 3;        //Number of spacial dimensions
	int fieldDim = 3;        //1 = scalar, 2,3=vector field
	FieldType fieldType = FIELD_MAGNETIC;  //What type of field (e.g. electric, magnetic)
	FieldClass fieldClass = FIELDCLASS_UNIFORM;  //What class of field? (e.g. uniform, gradient, dipole etc.)
	bool symmetry[3] = { false, false, false }; //Whether to mirror the field across an axis x,y,z
	bool useCubicInterpolation = false; //Whether to use cubic interpolation IF an interpolated field
};

////////////////////////////////////////////////////////////////
/// class SRKField
/// Abstract base class for electric or magnetic fields
/// Based on treatment in Geant4
///
/// Author: Matthew Bales
///////////////////////////////////////////////////////////////
class SRKField
{
public:

	SRKField();
	SRKField(SRKFieldSettings inpFS);
	virtual ~SRKField()
	{
	};

	//All SRKFields must implement addFieldValue
	virtual void addFieldValue(const double globalPoint[4], double fieldArray[9])=0; ///Field is an a
	template <class floatlike1, class floatlike2>
	inline void getFieldValue(const floatlike1 globalPoint[4], floatlike2 fieldArray[9])
	{
		fieldArray[0] = fieldArray[1] = fieldArray[2] = fieldArray[3] = fieldArray[4] = fieldArray[5] = fieldArray[6] = fieldArray[7] = fieldArray[8] = 0.;
		addFieldValue(globalPoint, fieldArray);
	}

	inline double getLength(){ return fs.extents[0];}
	inline double getWidth(){ return fs.extents[1];}
	inline double getHeight(){return fs.extents[2];}

	inline FieldType getFieldType(){ return fs.fieldType;}

protected:


	SRKFieldSettings fs;  //Field settings for the field
	int g4FieldX, g4FieldY, g4FieldZ;  /// index for fieldArray to fill  (0,1,2 = magnetic  3,4,5 = electric  6,7,8 = gravity
	double scalingFactorWUnits; //The fs.scaling factor multiplied units appropriate to the type and class

};

#endif
