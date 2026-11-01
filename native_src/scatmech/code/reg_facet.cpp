











#include "facet.h"
#include "torrspar.h"
#include "flake.h"


namespace SCATMECH {


    void Register(const Facet_BRDF_Model* x)
    {
        static bool Models_Registered = false;
        if (!Models_Registered) {
            Models_Registered=true;

            Register_Model(Facet_BRDF_Model);
            Register_Model(Shadowed_Facet_BRDF_Model);
            Register_Model(Subsurface_Facet_BRDF_Model);
        }
    }


}

