











#ifndef SCATMECH_SUBSURFACE_PARTICLE_H
#define SCATMECH_SUBSURFACE_PARTICLE_H

#include "local.h"
#include "sphrscat.h"
#include "filmtran.h"


namespace SCATMECH {


    class Subsurface_Particle_BRDF_Model : public Local_BRDF_Model
    {
        public:
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,depth);
            DECLARE_PARAMETER(Free_Space_Scatterer_Ptr,scatterer);
            DECLARE_PARAMETER(StackModel_Ptr,stack);



            DECLARE_PARAMETER(double,alpha);

            DECLARE_PARAMETER(double,beta);



        protected:
            virtual void setup();
            virtual JonesMatrix jonesDSC();
            Matrix Euler;
    };


}


#endif

