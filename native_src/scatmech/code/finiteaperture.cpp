











#include "scatmech.h"
#include "finiteaperture.h"
#include "focussedbeam.h"
#include "askuser.h"


namespace SCATMECH {





    Finite_Aperture_Instrument_BRDF_Model::
    Finite_Aperture_Instrument_BRDF_Model()
    {
        model_cs = xyxy;
    }

    MuellerMatrix
    Finite_Aperture_Instrument_BRDF_Model::
    mueller()
    {
        SETUP();

        if (lambda!=model->get_lambda()) error("lambda!=model.lambda");
        if (type!=model->get_type()) error("type!=model.type");
        if (substrate.index(lambda)!=model->get_substrate().index(lambda)) error("substrate!=model.substrate");


        if (alpha==0.) return model->Mueller(thetai,thetas,phis,rotation,xyxy);

        MuellerMatrix m=MuellerZero();

        Circle_Integral ci(integralmode);

        for (int i=0; i<ci.n(); ++i) {
            double theta = ci.theta(i);
            double r = ci.r(i);
            double w = ci.w(i);

            r = r*alpha*deg;

            Vector v = four_angles(thetas,phis,r,theta);

            double thetas_ = acos(v.z);
            double phis_ = atan2(v.y,v.x);

            m += model->Mueller(thetai,thetas_,phis_,rotation,xyxy)*w*cos(thetas_);
        }
        return (m/cos(thetas));
    }





    Vector
    Finite_Aperture_Instrument_BRDF_Model::
    four_angles(double theta,double phi,double alpha,double beta)
    {
        return Vector(
                   cos(beta)*sin(alpha)*sin(phi) +
                   cos(phi)*(cos(theta)*sin(alpha)*sin(beta) + cos(alpha)*sin(theta)),
                   -(cos(beta)*cos(phi)*sin(alpha)) +
                   sin(phi)*(cos(theta)*sin(alpha)*sin(beta) + cos(alpha)*sin(theta)),
                   cos(alpha)*cos(theta) - sin(alpha)*sin(beta)*sin(theta));
    }

    DEFINE_MODEL(Finite_Aperture_Instrument_BRDF_Model,Instrument_BRDF_Model,
                 "A BRDF model evaluated with a finite collection aperture.");

    DEFINE_PTRPARAMETER(Finite_Aperture_Instrument_BRDF_Model,BRDF_Model_Ptr,model,"BRDF Model to be integrated","Microroughness_BRDF_Model",0xFF);

    DEFINE_PARAMETER(Finite_Aperture_Instrument_BRDF_Model,double,alpha,"Opening half angle of detector [deg]","0",0xFF);
    DEFINE_PARAMETER(Finite_Aperture_Instrument_BRDF_Model,int,integralmode,"Order of Gauss-Zernike integral (1-7)","3",0xFF);


}



