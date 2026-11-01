











#include "zernike.h"
#include "zernikeexpansion.h"
#include "bobvlieg.h"

using namespace std;
using namespace SCATMECH;

namespace SCATMECH {


static double az(int k,double phi)
{
	if (k>0) return cos(k*phi);
	if (k<0) return -sin(k*phi);
	else  return 1./sqrt(2.0);
}

void ZernikeExpansion_BRDF_Model::setup()
{

	BRDF_Model::setup();


	model_cs = xyxy;


	ifstream_with_comments cfile(coefficientfile.c_str());


	if (!cfile) error("Cannot open coefficient file.");


	coeff.clear();


	while (!cfile.eof()) {
		short i,j,n,m,k,l,p;
		double c;


		string line = cfile.getstrline();


		if (line.size()!=0) {


			istringstream sline(line.c_str());


			sline >> i;
			while (isspace(sline.peek())||sline.peek()==',') sline.get();
			sline >> j;
			while (isspace(sline.peek())||sline.peek()==',') sline.get();
			sline >> n;
			while (isspace(sline.peek())||sline.peek()==',') sline.get();
			sline >> m;
			while (isspace(sline.peek())||sline.peek()==',') sline.get();
			sline >> k;
			while (isspace(sline.peek())||sline.peek()==',') sline.get();
			sline >> l;
			while (isspace(sline.peek())||sline.peek()==',') sline.get();
			sline >> p;
			while (isspace(sline.peek())||sline.peek()==',') sline.get();
			sline >> c;


			if (sline.fail()) error("Error reading a line in coefficient file");


			coeff.push_back(Coefficient(i,j,n,m,k,l,p,c));
		}
	}
	return;
}

MuellerMatrix ZernikeExpansion_BRDF_Model::mueller()
{


	SETUP();



	throw_transmission();
	throw_backward();


	MuellerMatrix result = MuellerZero();
	double brdf = 0;


	double rhoi = sqrt(2.)*sin(thetai/2.);
	double rhor = sqrt(2.)*sin(thetas/2.);
	double phii = pi-rotation;
	double phir = phis-rotation;


	for (list<Coefficient>::iterator q=coeff.begin();q!=coeff.end();++q) {
		short i = q->i;
		short j = q->j;
		short n = q->n;
		short m = q->m;
		short k = q->k;
		short l = q->l;
		short p = q->p;
		double coeff = q->coeff;


		double power = p==0?1.:p==1?lambda:p==2?sqr(lambda):p==3?sqr(lambda):pow(lambda,p);

		if (i==1&&j==1) {

			double prefact = 0.5/pi*sqrt((n+1.)*(m+1.)/((n==0||(n==m&&l==0))?4.:((n==m||l==0)?2.:1.)))*power;
			brdf += coeff * prefact * (Radial_Zernike_Polynomial(n,l,rhoi)*Radial_Zernike_Polynomial(m,l,rhor) +
						 Radial_Zernike_Polynomial(m,l,rhoi)*Radial_Zernike_Polynomial(n,l,rhor)) * cos(l*(phir-phii));
		} else if (i==j) {

			double prefact = 1./sqrt(2.+2.*(n==m?1:0)*(k==l?1:0))*sqrt((n+1.)*(m+1.))/pi*power;
			result[i-1][j-1] += coeff * prefact * (Radial_Zernike_Polynomial(n,k,rhoi)*Radial_Zernike_Polynomial(m,l,rhor)*az(k,phii)*az(l,phir) +
												 Radial_Zernike_Polynomial(m,l,rhoi)*Radial_Zernike_Polynomial(n,k,rhor)*az(k,phir)*az(l,phii));
		} else if (i<j) {

			double prefact = sqrt((n+1.)*(m+1.))/pi*power;
			result[i-1][j-1] += coeff * prefact * Radial_Zernike_Polynomial(n,k,rhoi)*Radial_Zernike_Polynomial(m,l,rhor)*az(k,phii)*az(l,phir);
			result[j-1][i-1] += coeff * prefact * Radial_Zernike_Polynomial(n,k,rhor)*Radial_Zernike_Polynomial(m,l,rhoi)*az(k,phir)*az(l,phii);
		} else throw("i>j in the coefficient file");
	}


	result[0][0] = 1.;

	result[2][0] = -result[2][0];
	result[2][1] = -result[2][1];
	result[3][2] = -result[3][2];


	result *= brdf*scale;


	return result;
}




DEFINE_MODEL(ZernikeExpansion_BRDF_Model,BRDF_Model,"Model using the model of Koenderink and van Doorn, extended by Germer to include the Mueller matrix");
DEFINE_PARAMETER(ZernikeExpansion_BRDF_Model,string,coefficientfile,"File containing coefficients","",0xFF);
DEFINE_PARAMETER(ZernikeExpansion_BRDF_Model,double,scale,"Scale factor","1",0xFF);


}
