#ifndef SRKGlobalField_HH
#define SRKGlobalField_HH 1

//  SRKGlobalField - handles the global ElectroMagnetic and Gravity field
//
//  There is a single SRKGlobalField object.
//
//  The field from each individual element is given by a
//  SRKField object. Any number of overlapping SRKField
//  objects can be added to the global field. Any element that
//  represents an element with an EMG field must add the appropriate
//  SRKField to the global GlobalField object.

#include "SRKField.h"
#include "SRKInterpolatedField.h"
#include "SRKDipoleField.h"
#include "SRKUniformField.h"
#include "SRKGradientField.h"

#include <vector>
#include <string>

const int MAX_SRK_NUM_FIELDS = 1000;

class SRKGlobalField: public SRKField
{

public:

	SRKGlobalField();
	SRKGlobalField(const SRKGlobalField&);

	virtual ~SRKGlobalField();

	SRKGlobalField& operator=(const SRKGlobalField&);

	void setupArray();

	/// GetFieldValue() returns the field value at a given point[].
	/// field is really field[6]: Bx,By,Bz,Ex,Ey,Ez,Gx,Gy,Gz
	/// point[] is in global coordinates: x,y,z,t.
	void addFieldValue(const double* point, double* outField) const;

	/// addField() adds the SRKField object for a single
	///  to the global field.
	void addField(SRKField* f)
	{
		theFieldList.push_back(f);
	}

	/// clear() removes all SRKField-s from the global object,
	/// and destroys them. Used before the geometry is completely
	/// re-created.
	void clear();

	/// updates all field tracking objects and clear()
	void updateField();

	///Creates/loads/adds all fields from fieldSettingsToLoad vector
	void constructFields();

	void setCurrentFieldSettingsToModify(int inp);

    inline void setFieldPath(std::string inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].fieldFilePath=inp;}
    inline void setFieldHistName(std::string inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].histName=inp;}
    inline void setFieldExtents(TVector3 inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].extents=inp;}
    inline void setFieldCenterPos(TVector3 inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].centerPos=inp;}
    inline void setFieldDirection(TVector3 inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].direction=inp;}
    inline void setFieldMoment(TVector3 inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].moment=inp;}
    inline void setFieldScalingValue(double inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].scalingValue=inp;}
    inline void setFieldFrequency(double inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].frequency=inp;}
    inline void setFieldSpaceDim(int inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].spaceDim=inp;}
    inline void setFieldFieldDim(int inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].fieldDim=inp;}
    inline void setFieldAngleX(double inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].angleX=inp;}
    inline void setFieldAngleY(double inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].angleY=inp;}
    inline void setFieldAngleZ(double inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].angleZ=inp;}
    inline void setFieldClass(int inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].fieldClass=static_cast<FieldClass>(inp);}
    inline void setFieldType(int inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].fieldType=static_cast<FieldType>(inp);}
    inline void setFieldSymmetryX(bool inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].symmetry[0]=inp;}
    inline void setFieldSymmetryY(bool inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].symmetry[1]=inp;}
    inline void setFieldSymmetryZ(bool inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].symmetry[2]=inp;}
    inline void setFieldUseCubicInterpolation(bool inp){if(currentFieldSettingsToModify >= 0) fieldSettingsToLoad[currentFieldSettingsToModify].useCubicInterpolation=inp;}
    inline void setFieldOffset(TVector3 inp){if(currentFieldSettingsToModify >= 0) {fieldSettingsToLoad[currentFieldSettingsToModify].offset=inp;fieldSettingsToLoad[currentFieldSettingsToModify].centerPos+=inp;}}


protected:

	static SRKGlobalField* object;

	int numFields;
	bool first;

	std::vector<SRKField*> theFieldList;

	SRKField** theFieldArray;

	std::vector<FieldSettings> fieldSettingsToLoad; //Settings to load magnetic or electric fields (only created during initialization);
	int currentFieldSettingsToModify;  //Current fieldSettingsToLoad to modify with macro commands...stupid macro command limitations

};

#endif
