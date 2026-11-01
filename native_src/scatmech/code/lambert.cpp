












#include "scatmech.h"
#include "lambert.h"

using namespace std;


namespace SCATMECH {


    MuellerMatrix
    Lambertian_BRDF_Model::
    mueller()
    {
        SETUP();

        return MuellerDepolarizer(reflectance->Get_Reflectance(lambda)/pi);
    }

    DEFINE_MODEL(Lambertian_BRDF_Model,BRDF_Model,
                 "A totally diffuse and depolarizing BRDF.");

    DEFINE_PTRPARAMETER(Lambertian_BRDF_Model,Reflectance_Ptr,reflectance,"Reflectance","Table_Reflectance",0xFF);


}


