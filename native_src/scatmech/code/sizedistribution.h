











#ifndef SCATMECH_SIZEDISTRIBUTION_H
#define SCATMECH_SIZEDISTRIBUTION_H

#include "inherit.h"

namespace SCATMECH {





	class Distribution : public Model {
	public:

		DECLARE_MODEL();

		virtual double pdf(double d) = 0;
	};

	typedef Model_Ptr<Distribution> Distribution_Ptr;
	void Register(const Distribution* x);




	class VolumeParticleSizeDistribution : public Model {
	public:

		virtual double volumedensity(double D) = 0;

		DECLARE_MODEL();
	};

	typedef Model_Ptr<VolumeParticleSizeDistribution> VolumeParticleSizeDistribution_Ptr;
	void Register(const VolumeParticleSizeDistribution* x);




	class Regular_VolumeParticleSizeDistribution : public Model {
	public:

		virtual double volumedensity(double D) {
			SETUP();
			return distribution->pdf(D)*numberdensity;
		};
		DECLARE_MODEL();
		DECLARE_PARAMETER(Distribution_Ptr, distribution);
		DECLARE_PARAMETER(double, numberdensity);
	};




	class SurfaceParticleSizeDistribution : public Model {
	public:

		virtual double surfacedensity(double d) = 0;

		DECLARE_MODEL();
	};

	typedef Model_Ptr<SurfaceParticleSizeDistribution> SurfaceParticleSizeDistribution_Ptr;
	void Register(const SurfaceParticleSizeDistribution* x);




	class Regular_SurfaceParticleSizeDistribution : public SurfaceParticleSizeDistribution {
	public:

		virtual double surfacedensity(double d) {
			SETUP();
			return distribution->pdf(d)*numberdensity;
		}

		DECLARE_MODEL();
		DECLARE_PARAMETER(Distribution_Ptr, distribution);
		DECLARE_PARAMETER(double, numberdensity);
	};




	class Log_Normal_Distribution : public Distribution {
	public:
		virtual double pdf(double d);

		DECLARE_MODEL();
		DECLARE_PARAMETER(double, sigma);
		DECLARE_PARAMETER(double, median);
	};




	class Normal_Distribution : public Distribution {
	public:
		virtual double pdf(double d);

		DECLARE_MODEL();
		DECLARE_PARAMETER(double, mean);
		DECLARE_PARAMETER(double, sigma);
	};




	class Weibull_Distribution : public Distribution {
	public:
		virtual double pdf(double d);

		DECLARE_MODEL();
		DECLARE_PARAMETER(double, exponent);
		DECLARE_PARAMETER(double, scale);
	};




	class Modified_Gamma_Distribution : public Distribution {
	public:
		virtual double pdf(double d);

		DECLARE_MODEL();
		DECLARE_PARAMETER(double, Lambda);
		DECLARE_PARAMETER(double, mu);
		DECLARE_PARAMETER(double, gamma);
	};




	class Bimodal_Distribution : public Distribution {
	public:
		virtual double pdf(double d);

		DECLARE_MODEL();
		DECLARE_PARAMETER(Distribution_Ptr, distributionA);
		DECLARE_PARAMETER(Distribution_Ptr, distributionB);
		DECLARE_PARAMETER(double, fractionA);
	};




	class CC1246E_SurfaceParticleSizeDistribution : public SurfaceParticleSizeDistribution {
	public:
		virtual double surfacedensity(double d);

		DECLARE_MODEL();
		DECLARE_PARAMETER(double, cleanliness);
		DECLARE_PARAMETER(double, slope);
	protected:
		void setup();
	private:
		double coverage(double slope, double level);
		double cov926;
		double cov;
	};

}

#endif
