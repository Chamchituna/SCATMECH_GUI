












#include <string>
#include "scatmech.h"
#include "urough.h"

using namespace std;


namespace SCATMECH {





    JonesMatrix
    Microroughness_BRDF_Model::
    jones()
    {

        SETUP();

        if (is_reflection()) {








            COMPLEX e = (COMPLEX)(substrate.epsilon(lambda));
            double _lambda = lambda;

            if (is_backward()) {
                e = 1./e;
                _lambda = _lambda/substrate.n(lambda);
            }


            COMPLEX pp = (e-1.)*
                         (
                             e*sin(thetai)*sin(thetas)-
                             sqrt(e-sqr(sin(thetas)))*sqrt(e-sqr(sin(thetai)))*cos(phis)
                         )/
                         (e*cos(thetai)+sqrt(e-sqr(sin(thetai))))/
                         (e*cos(thetas)+sqrt(e-sqr(sin(thetas))));


            COMPLEX ss = (e-1.)*cos(phis)/
                         (cos(thetai)+sqrt(e-sqr(sin(thetai))))/
                         (cos(thetas)+sqrt(e-sqr(sin(thetas))));


            COMPLEX ps = -(e-1.)*sqrt(e-sqr(sin(thetai)))*sin(phis)/
                         (e*cos(thetai)+sqrt(e-sqr(sin(thetai))))/
                         (cos(thetas)+sqrt(e-sqr(sin(thetas))));

            COMPLEX sp = -(e-1.)*sqrt(e-sqr(sin(thetas)))*sin(phis)/
                         (cos(thetai)+sqrt(e-sqr(sin(thetai))))/
                         (e*cos(thetas)+sqrt(e-sqr(sin(thetas))));


            double fx,fy;
            Bragg_Frequency(fx,fy);


            COMPLEX factor=sqrt(
                               16.*sqr(pi)/sqr(sqr(_lambda))*
                               cos(thetai)*cos(thetas)*
                               psd->psd(fx,fy)
                           );


            JonesMatrix result = JonesMatrix(pp,ss,ps,sp)*factor;
            if (is_forward()) {
                return result;
            } else {
                return JonesMatrix(result[0],result[1],-result[2],-result[3]);
            }

        } else {








            COMPLEX q0[2];
            COMPLEX q1[2];

            COMPLEX eps=(COMPLEX)(substrate.epsilon(lambda));
            double _lambda = lambda;

            if (is_backward()) {
                eps = 1./eps;
                _lambda = lambda/substrate.n(lambda);
            }

            COMPLEX index = sqrt(eps);
            COMPLEX k = sin(thetas)*index;
            double k0 = sin(thetai);

            q0[0]=sqrt(eps-k0*k0);
            q1[0]=sqrt(eps-k*k);
            q0[1]=sqrt(1.-k0*k0);
            q1[1]=sqrt(1.-k*k);

            double fx,fy;
            Bragg_Frequency(fx,fy);

            double PSD = psd->psd(fx,fy);

            JonesMatrix J;
            J.PP() =  index*(k*k0+cos(phis)*q0[0]*q1[1])/(q0[0]+eps*q0[1])/(q1[0]+eps*q1[1]);
            J.SS() =  cos(phis)/(q0[0]+q0[1])/(q1[0]+q1[1]);
            J.SP() =  index*q1[1]*sin(phis)/(q0[0]+q0[1])/(q1[0]+eps*q1[1]);
            J.PS() =  -q0[0]*sin(phis)/(q0[0]+eps*q0[1])/(q1[0]+q1[1]);

            COMPLEX common = 4.*pi*(1.-eps)*index*sqrt(q1[0]*q0[1]*PSD)/sqr(_lambda);

            if (is_forward()) {
                return J*common;
            } else {
                return JonesMatrix(J[0],J[1],-J[2],-J[3])*common;
            }
        }
    }

    DEFINE_MODEL(Microroughness_BRDF_Model,
                 Roughness_BRDF_Model,
                 "Scattering by a single rough surface in the smooth-surface limit (no films).");


}

