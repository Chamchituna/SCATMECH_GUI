












#ifndef SCATMECH_CROUGH_H
#define SCATMECH_CROUGH_H

#include "rough.h"
#include "allrough.h"


namespace SCATMECH {





    class Correlated_Roughness_BRDF_Model : public Roughness_BRDF_Model
    {
        public:

            DECLARE_MODEL();

            DECLARE_PARAMETER(double,thickness);

            DECLARE_PARAMETER(dielectric_function,film);

        protected:

            void setup();


            JonesMatrix jones();

        private:
            Correlated_Roughness_Stack_BRDF_Model model;

    };

}


#endif

