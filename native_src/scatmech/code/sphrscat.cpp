











#include "sphrscat.h"
#include "askuser.h"
#include "matrix3d.h"

using namespace std;

namespace SCATMECH {

	MuellerMatrix
	Free_Space_Scatterer::
	extinction(const Vector& v)
	{





		SETUP();

		JonesMatrix J = jones(v, v);

		MuellerMatrix result;
		result[0][0] = real(J.SS() + J.PP());
		result[0][1] = real(J.SS() - J.PP());
		result[0][2] = real(J.PS() + J.SP());
		result[0][3] = imag(J.SP() - J.PS());

		result[1][0] = result[0][1];
		result[1][1] = result[0][0];
		result[1][2] = real(J.PS() - J.SP());
		result[1][3] = -imag(J.PS() + J.SP());

		result[2][0] = result[0][2];
		result[2][1] = -result[1][2];
		result[2][2] = result[0][0];
		result[2][3] = imag(J.SS() - J.PP());

		result[3][0] = result[0][3];
		result[3][1] = -result[1][3];
		result[3][2] = -result[2][3];
		result[3][3] = result[0][0];

		double k = 2 * pi / lambda * medium.n(lambda);
		return (2 * pi / sqr(k)) * result;
	}

    DEFINE_VIRTUAL_MODEL(Free_Space_Scatterer,Model,
                         "Generalized free-space scatterer");

    DEFINE_PARAMETER(Free_Space_Scatterer,double,lambda,"Wavelength [um]","0.532",0xFF);
    DEFINE_PARAMETER(Free_Space_Scatterer,dielectric_function,medium,"Optical properties of the surrounding medium","(1,0)",0xFF);

    DEFINE_VIRTUAL_MODEL(SphericalScatterer,Free_Space_Scatterer,
                         "A spherically-symmetric free-space scatterer");

    DEFINE_PARAMETER(SphericalScatterer,double,radius,"Radius [um]","0.05",0xFF);
    DEFINE_PARAMETER(SphericalScatterer,dielectric_function,sphere,"Optical properties of the sphere","(1.59,0)",0xFF);


}



