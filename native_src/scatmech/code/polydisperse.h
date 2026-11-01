











#ifndef SCATMECH_POLYDISPERSE_H
#define SCATMECH_POLYDISPERSE_H

#include "bobvlieg.h"
#include "firstdiffuse.h"
#include <vector>
#include "sphrscat.h"
#include "oasphere.h"
#include "sphprt.h"

namespace SCATMECH {







	class Polydisperse_Sphere_BRDF_Model : public BRDF_Model {
	public:
		DECLARE_MODEL();


		DECLARE_PARAMETER(dielectric_function, particle);


		DECLARE_PARAMETER(StackModel_Ptr, stack);


		DECLARE_PARAMETER(SurfaceParticleSizeDistribution_Ptr,distribution);


		DECLARE_PARAMETER(double, Dstart);


		DECLARE_PARAMETER(double, Dend);


		DECLARE_PARAMETER(double, Dstep);



		DECLARE_PARAMETER(double, fractional_coverage);





		DECLARE_PARAMETER(double, antirainbow);

	protected:
		MuellerMatrix mueller();

	private:
		void setup();

		typedef Model_Ptr<Double_Interaction_BRDF_Model> Double_Interaction_BRDF_Model_Ptr;
		typedef std::list<Double_Interaction_BRDF_Model_Ptr> ModelList;
		typedef ModelList::iterator ModelListIterator;

		ModelList submodels;

		double distributioncoverage;
	};


}

#endif
