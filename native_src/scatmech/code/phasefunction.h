











#ifndef SCATMECH_PHASEFUNCTION_H
#define SCATMECH_PHASEFUNCTION_H

#include "inherit.h"
#include "sphrscat.h"
#include "sizedistribution.h"
#include "miescat.h"

namespace SCATMECH {




	class Phase_Function : public Model
	{
	public:



		virtual double f(double theta) = 0;

		DECLARE_MODEL();
	};

	typedef Model_Ptr<Phase_Function> Phase_Function_Ptr;
	void Register(const Phase_Function* x);






	class Polarized_Phase_Function : public Model {
	public:





		virtual MuellerMatrix f(const Vector& in, const Vector& out) = 0;



		virtual MuellerMatrix absorption_coefficient(const Vector& in) = 0;



		virtual double scattering_coefficient(const Vector& in) = 0;

		DECLARE_MODEL();
	};

	typedef Model_Ptr<Polarized_Phase_Function> Polarized_Phase_Function_Ptr;
	void Register(const Polarized_Phase_Function* x);









	class Henyey_Greenstein_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta);
		DECLARE_MODEL();
		DECLARE_PARAMETER(double, g);
	};





	class Double_Henyey_Greenstein_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta);
		DECLARE_MODEL();
		DECLARE_PARAMETER(double, g1);
		DECLARE_PARAMETER(double, g2);
		DECLARE_PARAMETER(double, c);
	};




	class Isotropic_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta) {
			return 0.25 / pi;
		}
		DECLARE_MODEL();
	};





	class Rayleigh_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta) {
			return 0.75*(1. + sqr(cos(theta)));
		}
		DECLARE_MODEL();
	};





	class Table_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta);
		DECLARE_MODEL();
		DECLARE_PARAMETER(Table, table);
	};





	class Reynolds_McCormick_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta);
		DECLARE_MODEL();
		DECLARE_PARAMETER(double, g);
		DECLARE_PARAMETER(double, alpha);
	protected:
		void setup();
	private:
		double K;
	};






	class Double_Reynolds_McCormick_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta);
		DECLARE_MODEL();
		DECLARE_PARAMETER(double, g);
		DECLARE_PARAMETER(double, alpha);
		DECLARE_PARAMETER(double, c);
	protected:
		void setup();
	private:
		double K;
	};




	class Legendre_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta);
		DECLARE_MODEL();
		DECLARE_PARAMETER(double, c0);
		DECLARE_PARAMETER(double, c1);
		DECLARE_PARAMETER(double, c2);
		DECLARE_PARAMETER(double, c3);
		DECLARE_PARAMETER(double, c4);
		DECLARE_PARAMETER(double, c5);
	};




	class Gaussian_Phase_Function : public Phase_Function
	{
	public:
		virtual double f(double theta);
		DECLARE_MODEL();
		DECLARE_PARAMETER(double, width);
		DECLARE_PARAMETER(double, background);
	};










	class Unpolarized_Phase_Function : public Polarized_Phase_Function {
	public:
		virtual MuellerMatrix f(const Vector& in, const Vector& out);
		virtual MuellerMatrix absorption_coefficient(const Vector& in);
		virtual double scattering_coefficient(const Vector& in);

		DECLARE_MODEL();
		DECLARE_PARAMETER(Phase_Function_Ptr, phase_function);
		DECLARE_PARAMETER(double, log_absorption);
		DECLARE_PARAMETER(double, log_scattering);
	};








	class CodeWhitney_Phase_Function : public Polarized_Phase_Function {
	public:
		virtual MuellerMatrix f(const Vector& in, const Vector& out);
		virtual MuellerMatrix absorption_coefficient(const Vector& in);
		virtual double scattering_coefficient(const Vector& in);

		DECLARE_MODEL();
		DECLARE_PARAMETER(Phase_Function_Ptr, phase_function);
		DECLARE_PARAMETER(double, log_absorption);
		DECLARE_PARAMETER(double, log_scattering);
		DECLARE_PARAMETER(double, plmax);
		DECLARE_PARAMETER(double, pcmax);
		DECLARE_PARAMETER(double, s);
	};









	class Constrained_Phase_Function : public Polarized_Phase_Function {
	public:
		virtual MuellerMatrix f(const Vector& in, const Vector& out);
		virtual MuellerMatrix absorption_coefficient(const Vector& in);
		virtual double scattering_coefficient(const Vector& in);

		DECLARE_MODEL();
		DECLARE_PARAMETER(Phase_Function_Ptr, phase_function);
		DECLARE_PARAMETER(double, log_absorption);
		DECLARE_PARAMETER(double, log_scattering);
		DECLARE_PARAMETER(double, m11_0);
		DECLARE_PARAMETER(double, m11_60);
		DECLARE_PARAMETER(double, m11_120);
		DECLARE_PARAMETER(double, m11_180);
		DECLARE_PARAMETER(double, m22_60);
		DECLARE_PARAMETER(double, m22_120);
		DECLARE_PARAMETER(double, m33_0);
		DECLARE_PARAMETER(double, m33_60);
		DECLARE_PARAMETER(double, m33_120);
		DECLARE_PARAMETER(double, m33_180);
		DECLARE_PARAMETER(double, m10_60);
		DECLARE_PARAMETER(double, m10_120);
		DECLARE_PARAMETER(double, m23_60);
		DECLARE_PARAMETER(double, m23_120);
	};





	class SphericalScatterer_Phase_Function : public Polarized_Phase_Function {
	public:
		virtual MuellerMatrix f(const Vector& in, const Vector& out);
		virtual MuellerMatrix absorption_coefficient(const Vector& in);
		virtual double scattering_coefficient(const Vector& in);
	protected:
		void setup();

	private:
		double csection;
		double absorption;
		double scattering;
		double k;

		DECLARE_MODEL();
		DECLARE_PARAMETER(SphericalScatterer_Ptr, scatterer);
		DECLARE_PARAMETER(double, volume_density);
	};





	class Distributed_Mie_Phase_Function : public Polarized_Phase_Function {
	public:
		virtual MuellerMatrix f(const Vector& in, const Vector& out);
		virtual MuellerMatrix absorption_coefficient(const Vector& in);
		virtual double scattering_coefficient(const Vector& in);
	protected:
		void setup();

		DECLARE_MODEL();
		DECLARE_PARAMETER(double, lambda);
		DECLARE_PARAMETER(dielectric_function, medium);
		DECLARE_PARAMETER(dielectric_function, sphere);
		DECLARE_PARAMETER(VolumeParticleSizeDistribution_Ptr, distribution);
		DECLARE_PARAMETER(double, integralStart);
		DECLARE_PARAMETER(double, integralEnd);
		DECLARE_PARAMETER(double, integralStep);

	private:
		typedef Model_Ptr<MieScatterer> MieScatterer_Ptr;
		typedef std::list<MieScatterer_Ptr> MieList;
		typedef MieList::iterator MieListIterator;
		typedef std::list<double> DensityList;
		typedef DensityList::iterator DensityListIterator;

		MieList mielist;
		DensityList densitylist;
		double norm;
		double scattering;
		double absorption;
	};

}
#endif
