











#ifndef ZERNIKEEXPANSION_H
#define ZERNIKEEXPANSION_H

#include "brdf.h"
#include <list>

namespace SCATMECH {








	class ZernikeExpansion_BRDF_Model : public BRDF_Model {
	public:


		DECLARE_MODEL();
		DECLARE_PARAMETER(std::string,coefficientfile);
		DECLARE_PARAMETER(double,scale);


		MuellerMatrix mueller();

	protected:


		void setup();


		struct Coefficient {
			Coefficient(short _i,short _j,short _n,short _m,short _k,short _l,short _p, double _coeff) : i(_i),j(_j),n(_n),m(_m),k(_k),l(_l),p(_p),coeff(_coeff) {}
			short i,j,n,m,k,l,p;
			double coeff;
		};


		std::list<Coefficient> coeff;
	};

}

#endif
