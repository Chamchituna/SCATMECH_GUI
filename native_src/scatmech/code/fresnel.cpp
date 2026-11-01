












#include "scatmech.h"
#include "fresnel.h"

using namespace std;


namespace SCATMECH {




    COMPLEX
    rs(dielectric_constant epsilon, double theta)
    {
        COMPLEX e=epsilon;
        return (cos(theta)-sqrt(e-sqr(sin(theta))))/
               (cos(theta)+sqrt(e-sqr(sin(theta))));
    }




    COMPLEX
    rp(dielectric_constant epsilon, double theta)
    {
        COMPLEX e=epsilon;
        return (e*cos(theta)-sqrt(e-sqr(sin(theta))))/
               (e*cos(theta)+sqrt(e-sqr(sin(theta))));
    }




    COMPLEX
    ts(dielectric_constant epsilon, double theta)
    {
        return 2.*cos(theta)/
               (cos(theta)+sqrt((COMPLEX)epsilon-sqr(sin(theta))));
    }




    COMPLEX
    tp(dielectric_constant epsilon, double theta)
    {
        return 2.*sqrt((COMPLEX)epsilon)*cos(theta)/
               ((COMPLEX)epsilon*cos(theta)+sqrt((COMPLEX)epsilon-sqr(sin(theta))));
    }




    double
    Ts(dielectric_constant epsilon, double theta)
    {
        return abs(sqrt((COMPLEX)epsilon-sqr(sin(theta)))/
                   cos(theta)*sqr(abs(ts(epsilon,theta))));
    }




    double
    Tp(dielectric_constant epsilon, double theta)
    {
        return abs(sqrt((COMPLEX)epsilon-sqr(sin(theta)))/
                   cos(theta)*sqr(abs(tp(epsilon,theta))));
    }









    COMPLEX
    sin_internal_angle(optical_constant n,COMPLEX thetai)
    {
        return sin(thetai)/(COMPLEX)n;
    }

    COMPLEX
    cos_internal_angle(optical_constant n,COMPLEX thetai)
    {
        return sqrt(COMPLEX(1,0)-sqr(sin(thetai)/(COMPLEX)(n)));
    }


}


