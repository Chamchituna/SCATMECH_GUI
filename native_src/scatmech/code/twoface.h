











#ifndef SCATMECH_TWOFACE_H
#define SCATMECH_TWOFACE_H

#include "scatmech.h"
#include "rough.h"
#include "roughnes.h"

namespace SCATMECH {







    class Two_Face_BRDF_Model: public Roughness_BRDF_Model
    {
        public:

            DECLARE_MODEL();


            DECLARE_PARAMETER(int,face);


            DECLARE_PARAMETER(double,thickness);


            DECLARE_PARAMETER(dielectric_function,film);

        protected:

            void setup();


            virtual JonesMatrix jones();

        private:
            Roughness_Stack_BRDF_Model model;

    };


}


#endif
