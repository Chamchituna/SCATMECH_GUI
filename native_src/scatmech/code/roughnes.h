












#ifndef SCATMECH_ROUGHNESS_H
#define SCATMECH_ROUGHNESS_H

#include "scatmech.h"
#include "urough.h"
#include "filmtran.h"


namespace SCATMECH {







    class Roughness_Stack_BRDF_Model :
        public Roughness_BRDF_Model
    {
        public:

            DECLARE_MODEL();


            DECLARE_PARAMETER(int,this_layer);


            DECLARE_PARAMETER(StackModel_Ptr,stack);

        protected:


            virtual JonesMatrix jones();

    };

    void Register(const Roughness_Stack_BRDF_Model* x);


}



#endif

