












#ifndef SCATMECH_UROUGH_H
#define SCATMECH_UROUGH_H

#include "scatmech.h"
#include "rough.h"


namespace SCATMECH {






    class Microroughness_BRDF_Model : public Roughness_BRDF_Model
    {
        public:
            DECLARE_MODEL();
        protected:

            JonesMatrix jones();
    };


}


#endif

