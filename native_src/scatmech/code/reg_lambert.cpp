











#include "lambert.h"
#include "diffuse.h"


namespace SCATMECH {

    void Register(const Lambertian_BRDF_Model* x)
    {
        static bool Models_Registered = false;
        if (!Models_Registered) {
            Models_Registered=true;

            Register_Model(Lambertian_BRDF_Model);
            Register_Model(Diffuse_Subsurface_BRDF_Model);
        }
    }


}

