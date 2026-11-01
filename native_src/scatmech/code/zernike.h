











#ifndef SCATMECH_ZERNIKE_H
#define SCATMECH_ZERNIKE_H

namespace SCATMECH {

	double Radial_Zernike_Polynomial(int n,int m,double rho);

	double Zernike_Polynomial(int n,int m,double x,double y);
	double Sin_Zernike_Polynomial(int n,int m,double x,double y);
	double Cos_Zernike_Polynomial(int n,int m,double x,double y);

};


#endif
