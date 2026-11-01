











#include "rough.h"
#include "urough.h"
#include "crough.h"
#include "roughnes.h"
#include "twoface.h"
#include "allrough.h"


namespace SCATMECH {


    void Register(const Roughness_BRDF_Model* x)
    {
        static bool Models_Registered = false;
        if (!Models_Registered) {
            Models_Registered=true;

            Register_Model(Roughness_BRDF_Model);
            Register_Model(Microroughness_BRDF_Model);
            Register_Model(Correlated_Roughness_BRDF_Model);
            Register((Roughness_Stack_BRDF_Model*)0);
            Register_Model(Two_Face_BRDF_Model);
        }
    }



}

