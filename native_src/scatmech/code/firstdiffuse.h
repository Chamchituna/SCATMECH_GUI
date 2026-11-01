











#ifndef SCATMECH_FIRSTDIFFUSE_H
#define SCATMECH_FIRSTDIFFUSE_H

#include "brdf.h"
#include "askuser.h"
#include "filmtran.h"
#include "phasefunction.h"


namespace SCATMECH {






	class First_Diffuse_BRDF_Model : public BRDF_Model
	{
	public:
		DECLARE_MODEL();
		DECLARE_PARAMETER(Polarized_Phase_Function_Ptr, phase_function);
		DECLARE_PARAMETER(StackModel_Ptr, stack);
		DECLARE_PARAMETER(double, depoll);
		DECLARE_PARAMETER(double, depolc);

	protected:
		virtual MuellerMatrix mueller();
	};

	void Register(const First_Diffuse_BRDF_Model*);

}

#endif

