











#include "polydisperse.h"
#include "miescat.h"

using namespace std;

namespace SCATMECH {

void
Polydisperse_Sphere_BRDF_Model::
setup()
{
    BRDF_Model::setup();


	submodels.clear();


	double sumdistribution = 0;
	double sumdistributionA = 0;


	for (double diameter = Dstart; diameter <= Dend; diameter *= 1. + Dstep) {


		submodels.push_back(new Double_Interaction_BRDF_Model);


		submodels.back()->set_lambda(lambda);
		submodels.back()->set_substrate(substrate);
		submodels.back()->set_type(type);



		double density = distribution->surfacedensity(diameter) * diameter * Dstep;
		submodels.back()->set_density(density);

		sumdistribution += density;
		sumdistributionA += density*pi*sqr(diameter / 2);

		double radius = diameter / 2;
		submodels.back()->set_stack(stack);
		submodels.back()->set_distance(radius);

		MieScatterer *mie = new MieScatterer;
		mie->set_medium(optical_constant(1, 0));
		mie->set_lambda(lambda);
		mie->set_radius(radius);



		optical_constant nsphere = particle.index(lambda);


		double absmin = antirainbow / diameter;


		double kmin = absmin / 4. / pi*lambda;


		if (nsphere.k < kmin) nsphere.k = kmin;
		mie->set_sphere(nsphere);

		submodels.back()->set_scatterer(mie);

		submodels.back()->set_alpha(0.);
		submodels.back()->set_beta(0.);
	}

	distributioncoverage = sumdistributionA;
}

MuellerMatrix
Polydisperse_Sphere_BRDF_Model::
mueller()
{
	SETUP();


	MuellerMatrix result = MuellerZero();


	for (ModelListIterator m = submodels.begin(); m != submodels.end(); ++m) {
		result += (*m)->Mueller(thetai, thetas, phis, rotation, model_cs);
	}

	if (fractional_coverage != 0) result *= fractional_coverage / distributioncoverage;

	return result;
}

DEFINE_MODEL(Polydisperse_Sphere_BRDF_Model,BRDF_Model,"Double interaction theory with a Mie scatterer and polydispersity.");
DEFINE_PTRPARAMETER(Polydisperse_Sphere_BRDF_Model,SurfaceParticleSizeDistribution_Ptr,distribution,"Size distribution","Regular_SurfaceParticleSizeDistribution",0xFF);
DEFINE_PTRPARAMETER(Polydisperse_Sphere_BRDF_Model,StackModel_Ptr,stack,"Films on substrate","No_StackModel",0xFF);
DEFINE_PARAMETER(Polydisperse_Sphere_BRDF_Model,dielectric_function,particle,"Optical properties of the particle","(1.5,0.0)",0xFF);
DEFINE_PARAMETER(Polydisperse_Sphere_BRDF_Model,double,Dstart,"Starting diameter [um]","1",0xFF);
DEFINE_PARAMETER(Polydisperse_Sphere_BRDF_Model,double,Dend,"Ending diameter [um]","1000",0xFF);
DEFINE_PARAMETER(Polydisperse_Sphere_BRDF_Model,double,Dstep,"Fractional diameter step","0.01",0xFF);
DEFINE_PARAMETER(Polydisperse_Sphere_BRDF_Model, double, fractional_coverage, "Fractional area coverage (0 = use distribution's coverage)", "0", 0xFF);
DEFINE_PARAMETER(Polydisperse_Sphere_BRDF_Model, double, antirainbow, "Anti-rainbow factor", "0", 0xFF);

}
