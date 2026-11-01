











#include "sphprt.h"
#include "sphdfct.h"
#include "onelayer.h"
#include "bobvlieg.h"
#include "subsphere.h"
#include "axipart.h"
#include "raystack.h"
#include "subbobvlieg.h"

namespace SCATMECH {


    void Register(const Local_BRDF_Model* x)
    {
        static bool Models_Registered = false;
        if (!Models_Registered) {
            Models_Registered=true;

            Register_Model(Local_BRDF_Model);
            Register_Model(Rayleigh_Defect_BRDF_Model);
            Register_Model(Rayleigh_Stack_BRDF_Model);
            Register_Model(OneLayer_BRDF_Model);
            Register_Model(Double_Interaction_BRDF_Model);
            Register_Model(Subsurface_Particle_BRDF_Model);
            Register_Model(Bobbert_Vlieger_BRDF_Model);
            Register_Model(Axisymmetric_Particle_BRDF_Model);
            Register_Model(Subsurface_Bobbert_Vlieger_BRDF_Model);
            Register_Model(Subsurface_Axisymmetric_Particle_BRDF_Model);
            Register((Axisymmetric_Shape*)0);
        }
    }


}

