












#include "scatmech.h"
#include "diffuse.h"
#include "fresnel.h"

using namespace std;


namespace SCATMECH {


    static
    COMPLEX
    ArcSin(const COMPLEX& a)
    {
        COMPLEX temp = sqrt(1. - sqr(a)) + COMPLEX(0,1)*a;
        return COMPLEX(arg(temp),-log(abs(temp)));
    }

    Diffuse_Subsurface_BRDF_Model::
    Diffuse_Subsurface_BRDF_Model()
    {
    }

    MuellerMatrix Diffuse_Subsurface_BRDF_Model::
    mueller()
    {
        SETUP();

        throw_transmission();
        throw_backward();

        optical_constant vacuum(1.,0.);

        double s = abs(sqrt((COMPLEX)substrate.epsilon(lambda)-sqr(sin(thetai)))/cos(thetai)
                       * sqrt((COMPLEX)substrate.epsilon(lambda)-sqr(sin(thetas)))/cos(thetas));


        s *= feedback/sqr(substrate.n(lambda));

        MuellerMatrix m = Lambertian_BRDF_Model::mueller();

        return s*(((MuellerMatrix)(stack->t12(thetas,lambda,vacuum,substrate))*m)*
                  (MuellerMatrix)(stack->t12(thetai,lambda,vacuum,substrate)));
    }

    void
    Diffuse_Subsurface_BRDF_Model::
    setup()
    {
        Lambertian_BRDF_Model::setup();



        substrate.force_nonabsorbing();

        double rho=0;
        double normalization=0;
        for (double theta=0; theta<=pi/2; theta+=deg/5.) {

            COMPLEX ti = ArcSin((COMPLEX)(sin(theta)*substrate.n(lambda)));


            JonesMatrix r = stack->r21(ti,lambda,substrate,vacuum);
            double rr=0.5*(norm(r[0])+norm(r[1]));


            double dOmega = sin(2.*theta);
            rho += rr*dOmega;
            normalization += dOmega;
        }
        rho = rho/normalization;



        double R = reflectance->Get_Reflectance(lambda);
        feedback = 1./(1.-rho*R);
    }

    DEFINE_MODEL(Diffuse_Subsurface_BRDF_Model,Lambertian_BRDF_Model,
                 "Totally diffuse scattering under a smooth surface.");

    DEFINE_PTRPARAMETER(Diffuse_Subsurface_BRDF_Model,StackModel_Ptr,stack,"Film stack on substrate","No_StackModel",0xFF);


}




