#ifndef SRKField_h
#define SRKField_h 1

#include "SRKROOTField.h"

#include <string>

const static double inv_c_light=1./2.99792458e+8;

enum FieldType{FIELD_MAGNETIC,FIELD_ELECTRIC,FIELD_GRAVITY};
enum FieldClass{FIELDCLASS_UNIFORM,FIELDCLASS_GRADIENT,FIELDCLASS_DIPOLE,FIELDCLASS_OSCILLATING,FIELDCLASS_INTERPOLATION};

//Simple settings class
class FieldSettings
{
public:
	FieldSettings()
	{
		scalingValue = 0;
		fieldType = FIELD_MAGNETIC;
		fieldClass = FIELDCLASS_UNIFORM;
		spaceDim = fieldDim = 3;
		angleX = angleY = angleZ = 0;
		symmetry[0] = symmetry[1] = symmetry[2] = false;
		histName = "field";
		extents = TVector3(100., 100. , 100. );
		fieldFilePath = "";
		useCubicInterpolation = false;
		frequency = 0;
	}
	;

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
		return out;
	}
	std::string fieldFilePath;
	std::string histName;

	TVector3 centerPos;  //Center position of field extents. Base for local unit system local unit system for field. Not used by interpolated field.[mm]
	TVector3 extents;    //Extents of the field in both directions from the center position (later to be changed)
	TVector3 offset;     //Offset of local unit system
	TVector3 direction;  //Direction of field (changed to unit vector)
	double frequency; 		  //rad/s
	TVector3 moment;	  //

	double scalingValue; //[tesla][V/m][m/s2] Zero scaling value will be the signal that no field should be loaded;
	int spaceDim;
	int fieldDim;
	double angleX;
	double angleY;
	double angleZ;
	FieldType fieldType;
	FieldClass fieldClass;
	bool symmetry[3];
	bool useCubicInterpolation;
};

class SRKField
{
public:

	//  Constructors
	SRKField();
	SRKField(FieldSettings inpFS);

	//  Destructor.
	virtual ~SRKField();

	inline double getLength(){ return fs.extents[0];}


	inline double getWidth(){ return fs.extents[1];}

	inline double getHeight(){return fs.extents[2];}

	inline bool DoesFieldChangeEnergy() const{ return true;}

	//All UCCNTFields must implement addFieldValue  This is for extra speed
	virtual void addFieldValue(const double Point[4], double Bfield[]) const=0;

	inline void GetFieldValue(const double Point[4], double Bfield[]) const
	{
		Bfield[0] = Bfield[1] = Bfield[2] = Bfield[3] = Bfield[4] = Bfield[5] = Bfield[6] = Bfield[7] = Bfield[8] = 0.;
		addFieldValue(Point, Bfield);
	}

	inline FieldType getFieldType(){ return fs.fieldType;}

protected:

	FieldSettings fs;
	int g4FieldX, g4FieldY, g4FieldZ;
	double scalingFactorWUnits;

};

#endif
