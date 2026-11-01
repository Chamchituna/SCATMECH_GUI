












#ifndef SCATMECH_DIFFUSE_H
#define SCATMECH_DIFFUSE_H

#include "filmtran.h"
#include "lambert.h"


namespace SCATMECH {

    class Diffuse_Subsurface_BRDF_Model: public Lambertian_BRDF_Model
    {
        public:

            Diffuse_Subsurface_BRDF_Model();

        public:

            DECLARE_MODEL();

            DECLARE_PARAMETER(StackModel_Ptr,stack);

        protected:
            MuellerMatrix mueller();
            void setup();

        private:
            double feedback;
    };


}




#endif
