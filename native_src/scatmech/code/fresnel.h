











#ifndef FRESNEL_H
#define FRESNEL_H

#include "optconst.h"


namespace SCATMECH {




    COMPLEX rs( dielectric_constant epsilon, double theta );
    COMPLEX rp( dielectric_constant epsilon, double theta );




    inline double Rs( dielectric_constant epsilon, double theta )
    {
        using namespace std;
        COMPLEX _rs=rs(epsilon,theta);
        return abs(_rs*conj(_rs));
    }

    inline double Rp( dielectric_constant epsilon, double theta )
    {
        using namespace std;
        COMPLEX _rp=rp(epsilon,theta);
        return abs(_rp*conj(_rp));
    }




    COMPLEX ts( dielectric_constant epsilon, double theta );
    COMPLEX tp( dielectric_constant epsilon, double theta );




    double Ts( dielectric_constant epsilon, double thetas );
    double Tp( dielectric_constant epsilon, double thetas );




    COMPLEX sin_internal_angle( optical_constant n, COMPLEX thetai );
    COMPLEX cos_internal_angle( optical_constant n, COMPLEX thetai );

}


#endif

