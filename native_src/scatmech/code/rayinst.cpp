











#include "rayinst.h"
#include "filmtran.h"
#include "askuser.h"

using namespace std;


namespace SCATMECH {


    static double subtend(double theta1,double phi1,double theta2, double phi2)
    {


        double cosreturnvalue = cos(theta1)*cos(theta2) +
                                cos(phi1)*cos(phi2)*sin(theta1)*sin(theta2) +
                                sin(phi1)*sin(phi2)*sin(theta1)*sin(theta2);
        if (cosreturnvalue>1) return 0.;
        if (cosreturnvalue<-1) return pi;
        else return acos(cosreturnvalue);
    }




    MuellerMatrix
    Rayleigh_Instrument_BRDF_Model::
    mueller()
    {

        SETUP();

        throw_backward();
        throw_transmission();

        if (is_reflection()) {







            double n_minus_1 = air.n(lambda)-1.;


            complex<double> rss = stack->rs12(thetas,lambda,vacuum,substrate);
            complex<double> rsi = stack->rs12(thetai,lambda,vacuum,substrate);
            complex<double> rps = stack->rp12(thetas,lambda,vacuum,substrate);
            complex<double> rpi = stack->rp12(thetai,lambda,vacuum,substrate);


            double subtendi = subtend(-thetai,0.,thetas,phis);
            double sini = fabs(sin(subtendi));
            double subtendr = subtend(thetai,0.,thetas,phis);
            double sinr = fabs(sin(subtendr));


            JonesMatrix term1(
                (complex<double>)(sin(thetai)*sin(thetas)-cos(thetas)*cos(thetai)*cos(phis)),
                (complex<double>)cos(phis),
                (complex<double>)(-cos(thetai)*sin(phis)),
                (complex<double>)(-cos(thetas)*sin(phis)));

            JonesMatrix term2(
                (-rpi*sin(thetai)*sin(thetas)-rpi*cos(thetas)*cos(thetai)*cos(phis)),
                rsi*cos(phis),
                -rpi*cos(thetai)*sin(phis),
                -rsi*cos(thetas)*sin(phis));

            JonesMatrix term3(
                (-rps*sin(thetai)*sin(thetas)-rps*cos(thetas)*cos(thetai)*cos(phis)),
                -rss*cos(phis),
                rss*cos(thetai)*sin(phis),
                -rps*cos(thetas)*sin(phis));

            JonesMatrix term4(
                (rps*rpi*sin(thetai)*sin(thetas)-rps*rpi*cos(thetas)*cos(thetai)*cos(phis)),
                -rss*rsi*cos(phis),
                rss*rpi*cos(thetai)*sin(phis),
                -rps*rsi*cos(thetas)*sin(phis));


            MuellerMatrix sterm1 = term1;
            MuellerMatrix sterm2 = term2;
            MuellerMatrix sterm3 = term3;
            MuellerMatrix sterm4 = term4;


            double i0=4.*sqr(pi)/pow(lambda,4.)*sqr(n_minus_1)/number_density;


            double factor = i0/2./cos(thetas)*field_of_view;
            MuellerMatrix result=((sterm1/sini)+(sterm2/sinr)+(sterm3/sinr)+(sterm4/sini))*factor;
            return result;
        }
        return MuellerZero();
    }

    DEFINE_MODEL(Rayleigh_Instrument_BRDF_Model,Instrument_BRDF_Model,
                 "Effective BRDF due to Rayleigh scatter by the air surrounding a smooth sample.");

    DEFINE_PTRPARAMETER(Rayleigh_Instrument_BRDF_Model,StackModel_Ptr,stack,"Film stack on substrate","No_StackModel",0xFF);

    DEFINE_PARAMETER(Rayleigh_Instrument_BRDF_Model,double,field_of_view,"Field of view of detector [um]","1000",0xFF);

    DEFINE_PARAMETER(Rayleigh_Instrument_BRDF_Model,dielectric_function,air,"Index of refraction of air","1.0002748",0xFF);

    DEFINE_PARAMETER(Rayleigh_Instrument_BRDF_Model,double,number_density,"Number density of air [um^-3]","2.51E7",0xFF);



}

