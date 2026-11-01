












#ifndef SCATMECH_SPHDFCT_H
#define SCATMECH_SPHDFCT_H

#include "scatmech.h"
#include "filmtran.h"
#include "local.h"


namespace SCATMECH {







    class Rayleigh_Defect_BRDF_Model :
        public Local_BRDF_Model
    {
        public:

            DECLARE_MODEL();


            DECLARE_PARAMETER(double,radius);


            DECLARE_PARAMETER(double,distance);


            DECLARE_PARAMETER(dielectric_function,defect);

            DECLARE_PARAMETER(StackModel_Ptr,stack);

        protected:



            void setup();






            double kd;


            double kr;

        protected:

            virtual JonesMatrix jonesDSC();

    };


}



#endif
