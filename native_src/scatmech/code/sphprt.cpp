












#include "scatmech.h"
#include "sphprt.h"
#include "vector3d.h"
#include "matrix3d.h"
#include "fresnel.h"
#include "askuser.h"
#include "rayscat.h"

using namespace std;


namespace SCATMECH {

    static
    COMPLEX
    ArcSin(const COMPLEX& a)
    {
        COMPLEX temp = sqrt(1. - sqr(a)) + COMPLEX(0,1)*a;
        return COMPLEX(arg(temp),-log(abs(temp)));
    }


    JonesMatrix
    Double_Interaction_BRDF_Model::
    jonesDSC()
    {
        SETUP();

        if (scatterer->get_medium().index(lambda)!=optical_constant(1,0)) error("scaterer.medium != 1");
        if (scatterer->get_lambda()!=lambda) error("scaterer.lambda != lambda");

        Matrix _Euler = Euler*Matrix(cos(rotation),sin(rotation),0,
                                     -sin(rotation),cos(rotation),0,
                                     0,0,1);

        if (fabs(sin(thetai)-sin(thetas)*cos(phis))<1E-5) thetai+=2.5E-5;
        if (fabs(sin(thetai)+sin(thetas)*cos(phis))<1E-5) thetai+=2.5E-5;

        if (is_forward()) {
            if (is_reflection()) {











                double k = 2.*pi/lambda;
                double kd = k*distance;

                if (fabs(thetai)<1E-5) thetai=1E-5;
                if (fabs(thetas)<1E-5) thetas=1E-5;


                Vector inK(sin(thetai),0,-cos(thetai));
                Vector outK(sin(thetas)*cos(phis),sin(thetas)*sin(phis),cos(thetas));


                Vector normal(0.,0.,1.);


                Vector inKr(inK.x,inK.y,-inK.z);
                Vector outKr(outK.x,outK.y,-outK.z);



                Vector inS =  perpto(normal,inK);
                Vector inP =  perpto(inK,inS);
                Vector outS = perpto(normal,outK);
                Vector outP = perpto(outK,outS);


                Vector inSr =  inS;
                Vector inPr =  perpto(inKr,inSr);
                Vector outSr = outS;
                Vector outPr = perpto(outKr,outSr);




                Vector b1in = perpto(inK,outK),   b1out=b1in;
                Vector b2in = perpto(inKr,outK),  b2out=b2in;
                Vector b3in = perpto(inK,outKr),  b3out=b3in;
                Vector b4in = perpto(inKr,outKr), b4out=b4in;

                Vector a1in = cross(b1in,inK),   a1out = cross(b1out,outK);
                Vector a2in = cross(b2in,inKr),  a2out = cross(b2out,outK);
                Vector a3in = cross(b3in,inK),   a3out = cross(b3out,outKr);
                Vector a4in = cross(b4in,inKr),  a4out = cross(b4out,outKr);


                JonesMatrix ri = stack->r12(thetai,lambda,vacuum,substrate);
                JonesMatrix rs = stack->r12(thetas,lambda,vacuum,substrate);

                JonesMatrix S1 = scatterer->jones(_Euler*inK, _Euler*outK);
                JonesMatrix S2 = scatterer->jones(_Euler*inKr,_Euler*outK);
                JonesMatrix S3 = scatterer->jones(_Euler*inK, _Euler*outKr);
                JonesMatrix S4 = scatterer->jones(_Euler*inKr,_Euler*outKr);


                COMPLEX I(0,1);
                COMPLEX alpha = exp(2.*I*kd*cos(thetai));
                COMPLEX beta  = exp(2.*I*kd*cos(thetas));

                JonesMatrix RotIn1  = GetJonesRotator(a1in,b1in,inS,inP);
                JonesMatrix RotIn2  = GetJonesRotator(a2in,b2in,inSr,inPr);
                JonesMatrix RotIn3  = GetJonesRotator(a3in,b3in,inS,inP);
                JonesMatrix RotIn4  = GetJonesRotator(a4in,b4in,inSr,inPr);

                JonesMatrix RotOut1 = GetJonesRotator(a1out,b1out,outS,outP).transpose();
                JonesMatrix RotOut2 = GetJonesRotator(a2out,b2out,outS,outP).transpose();
                JonesMatrix RotOut3 = GetJonesRotator(a3out,b3out,outSr,outPr).transpose();
                JonesMatrix RotOut4 = GetJonesRotator(a4out,b4out,outSr,outPr).transpose();

                S1 =    RotOut1*S1*RotIn1;
                S2 =    RotOut2*S2*RotIn2*alpha*ri;
                S3 = rs*RotOut3*S3*RotIn3*beta;
                S4 = rs*RotOut4*S4*RotIn4*alpha*beta*ri;

                COMPLEX common=1./k;

                return (S1+S2+S3+S4)*common;

            } else {
                double n = substrate.n(lambda);

                double k=2.*pi/lambda;

                Vector normal(0,0,1);


                Vector inK(sin(thetai),0.,-cos(thetai));

                Vector inKr(inK.x,inK.y,-inK.z);


                Vector inS = perpto(inK,normal);
                Vector inP = cross(inK,inS);
                Vector inSr = inS;
                Vector inPr = cross(inKr,inSr);


                double sin_thetas_outside = sin(thetas)*n;


                if (fabs(sin_thetas_outside)>1.) return JonesZero();

                double thetas_outside = asin(sin_thetas_outside);
                double cos_thetas_outside = sqrt(1.-sqr(sin_thetas_outside));


                Vector outK(sin_thetas_outside*cos(phis),sin_thetas_outside*sin(phis),-cos_thetas_outside);
                Vector outS = perpto(outK,normal);
                Vector outP = cross(outK,outS);


                JonesMatrix scatter_direct;
                {

                    Vector perp = perpto(inK,outK);
                    Vector pari = cross(perp,inK);
                    Vector pars = cross(perp,outK);


                    JonesMatrix matrixin = GetJonesRotator(pari,perp,inS,inP);
                    JonesMatrix matrixout = GetJonesRotator(outS,outP,pars,perp);


                    JonesMatrix scatter = scatterer->jones(_Euler*inK,_Euler*outK);
                    scatter_direct = matrixout*scatter*matrixin;
                }


                JonesMatrix scatter_indirect;
                {

                    Vector perp = perpto(inKr,outK);
                    Vector pari = cross(perp,inKr);
                    Vector pars = cross(perp,outK);


                    JonesMatrix matrixin = GetJonesRotator(pari,perp,inSr,inPr);
                    JonesMatrix matrixout = GetJonesRotator(outS,outP,pars,perp);


                    JonesMatrix scatter = scatterer->jones(_Euler*inKr,_Euler*outK);


                    COMPLEX phase = exp(COMPLEX(0,2)*cos(thetai)*k*distance);
                    JonesMatrix r = stack->r12(thetai,lambda,vacuum,substrate);

                    scatter_indirect = phase*matrixout*scatter*matrixin*r;
                }


                JonesMatrix t = stack->t12(thetas_outside,lambda,vacuum,substrate)*
                                sqrt(n*cos(thetas)/cos_thetas_outside);


                JonesMatrix scatter = t*(scatter_direct+scatter_indirect)/k;


                double jacobian = sqr(n)*cos(thetas)/cos_thetas_outside;


                COMPLEX phase = exp(COMPLEX(0, 1) * k * distance * (cos_thetas_outside - n * cos(thetas)));


                return scatter*sqrt(jacobian)*phase;
            }
        } else {

            double n = substrate.n(lambda);
            double k=2.*pi/lambda;

            if (is_reflection()) {


                double sin_thetai_at_part = sin(thetai)*n;


                if (sin_thetai_at_part>1) return JonesZero();

                double cos_thetai_at_part = sqrt(1.-sqr(sin_thetai_at_part));


                Vector inKi(sin_thetai_at_part,0.,cos_thetai_at_part);


                Vector normal(0.,0.,1.);
                Vector inSi =  perpto(inKi,normal);
                Vector inPi =  cross(inKi,inSi);


                double sin_thetas_at_part = sin(thetas)*n;


                if (sin_thetas_at_part>1) return JonesZero();

                double cos_thetas_at_part = sqrt(1.-sqr(sin_thetas_at_part));


                Vector outKi(sin_thetas_at_part*cos(phis),sin_thetas_at_part*sin(phis),-cos_thetas_at_part);


                Vector outSi = perpto(outKi,normal);
                Vector outPi = cross(outKi,outSi);


                JonesMatrix ti = stack->t21i(thetai,lambda,substrate,vacuum)/(COMPLEX)sqrt(n);
                JonesMatrix ts = stack->t21i(thetas,lambda,substrate,vacuum)*(COMPLEX)sqrt(cos_thetas_at_part/cos(thetas)/n);


                Vector perp = perpto(inKi,outKi);
                Vector pari = cross(perp,inKi);
                Vector pars = cross(perp,outKi);


                JonesMatrix matrixin = GetJonesRotator(pari,perp,inSi,inPi);
                JonesMatrix matrixout = GetJonesRotator(outSi,outPi,pars,perp);


                JonesMatrix scatter = scatterer->jones(_Euler*inKi,_Euler*outKi);
                scatter=matrixout*scatter*matrixin;


                double jacobian = cos(thetas)/cos_thetas_at_part*sqr(n);


                scatter = ts*scatter*ti;

                COMPLEX common=jacobian/sqr(k);


                COMPLEX phasei = exp(COMPLEX(0, 1) * k * distance * (cos_thetai_at_part - n * cos(thetai)));
                COMPLEX phases = exp(COMPLEX(0, 1) * k * distance * (cos_thetas_at_part - n * cos(thetas)));


                return scatter*sqrt(common)*phasei*phases;

            } else {


                double sin_thetai_at_part = sin(thetai)*n;


                if (sin_thetai_at_part>1) return JonesZero();

                double cos_thetai_at_part = sqrt(1.-sqr(sin_thetai_at_part));


                Vector inK(sin_thetai_at_part,0.,cos_thetai_at_part);


                Vector normal(0.,0.,1.);
                Vector inS =  perpto(inK,normal);
                Vector inP =  cross(inK,inS);

                if (abs(thetas)<1E-8) thetas=1E-8;


                double sin_thetas = sin(thetas);
                double cos_thetas = sqrt(1.-sqr(sin_thetas));


                Vector outKi(sin_thetas*cos(phis),sin_thetas*sin(phis),cos_thetas);
                Vector outSi = perpto(outKi,normal);
                Vector outPi = cross(outKi,outSi);


                Vector outKr(sin_thetas*cos(phis),sin_thetas*sin(phis),-cos_thetas);
                Vector outSr = outSi;
                Vector outPr = cross(outKr,outSr);


                JonesMatrix scatter_direct;
                {

                    Vector perp = perpto(inK,outKi);
                    Vector pari = cross(perp,inK);
                    Vector pars = cross(perp,outKi);


                    JonesMatrix matrixin = GetJonesRotator(pari,perp,inS,inP);
                    JonesMatrix matrixout = GetJonesRotator(outSi,outPi,pars,perp);


                    JonesMatrix scatter = scatterer->jones(_Euler*inK,_Euler*outKi);


                    scatter_direct=matrixout*scatter*matrixin;
                }


                JonesMatrix scatter_indirect;
                {

                    Vector perp = perpto(inK,outKr);
                    Vector pari = cross(perp,inK);
                    Vector pars = cross(perp,outKr);


                    JonesMatrix r = stack->r12(thetas,lambda,vacuum,substrate);


                    JonesMatrix matrixin = GetJonesRotator(pari,perp,inS,inP);
                    JonesMatrix matrixout = GetJonesRotator(outSr,outPr,pars,perp);


                    JonesMatrix scatter = scatterer->jones(_Euler*inK,_Euler*outKr);


                    COMPLEX phase = exp(COMPLEX(0,2)*cos_thetas*k*distance);

                    scatter_indirect=r*matrixout*scatter*matrixin*phase;
                }


                JonesMatrix ti = stack->t21i(thetai,lambda,substrate,vacuum)/(COMPLEX)sqrt(n);

                COMPLEX phasei = exp(COMPLEX(0, 1) * k * distance * (cos_thetai_at_part - n * cos(thetai)));


                JonesMatrix scatter = (scatter_direct+scatter_indirect)*ti/k*phasei;

                return scatter;

            }
        }
    }




    void
    Double_Interaction_BRDF_Model::
    setup()
    {

        Local_BRDF_Model::setup();
        if (scatterer->get_lambda()!=lambda) scatterer->set_lambda(lambda);
        if (scatterer->get_medium().n(lambda)!=1) scatterer->set_medium(vacuum);
        Euler = Matrix(cos(alpha)*cos(beta),  sin(alpha), -cos(alpha)*sin(beta),
                       -cos(beta)*sin(alpha), cos(alpha), sin(alpha)*sin(beta),
                       sin(beta),                 0,            cos(beta));
    }

    DEFINE_MODEL(Double_Interaction_BRDF_Model,Local_BRDF_Model,
                 "The double-interaction model for a spherical particle above a surface.");

    DEFINE_PTRPARAMETER(Double_Interaction_BRDF_Model,StackModel_Ptr,stack,"Film stack on substrate","No_StackModel",0xFF);
    DEFINE_PARAMETER(Double_Interaction_BRDF_Model,double,distance,"Distance from center to surface [um]","0.05",0xFF);
    DEFINE_PTRPARAMETER(Double_Interaction_BRDF_Model,Free_Space_Scatterer_Ptr,scatterer,"The scattering function","MieScatterer",0xFF);
    DEFINE_PARAMETER(Double_Interaction_BRDF_Model,double,alpha,"Rotation of particle about z-axis [rad]","0",0xFF);
    DEFINE_PARAMETER(Double_Interaction_BRDF_Model,double,beta,"Rotation of particle about y-axis [rad]","0",0xFF);


}


