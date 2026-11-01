












#ifndef SCATMECH_DIPARTICLE_H
#define SCATMECH_DIPARTICLE_H

#include "scatmech.h"
#include "filmtran.h"
#include "sphrscat.h"
#include "local.h"
#include "matrix3d.h"


namespace SCATMECH {













    class Double_Interaction_BRDF_Model :
        public Local_BRDF_Model
    {
        public:






        public:
            DECLARE_MODEL();



            DECLARE_PARAMETER(double,distance);


            DECLARE_PARAMETER(Free_Space_Scatterer_Ptr,scatterer);


            DECLARE_PARAMETER(StackModel_Ptr,stack);







            DECLARE_PARAMETER(double,alpha);
            DECLARE_PARAMETER(double,beta);

        protected:






            void setup();


            virtual JonesMatrix jonesDSC();


            Matrix Euler;
    };


}


#endif
