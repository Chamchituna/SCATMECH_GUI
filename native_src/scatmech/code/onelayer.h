











#ifndef SCATMECH_ONELAYER_H
#define SCATMECH_ONELAYER_H

#include "scatmech.h"
#include "local.h"


namespace SCATMECH {








    class OneLayer_BRDF_Model : public Local_BRDF_Model
    {
        public:

            DECLARE_MODEL();


            DECLARE_PARAMETER(double,radius);


            DECLARE_PARAMETER(dielectric_function,defect);


            DECLARE_PARAMETER(dielectric_function,film);


            DECLARE_PARAMETER(double,tau);



            DECLARE_PARAMETER(double,depth);

        protected:

            JonesMatrix jonesDSC();

    };


}


#endif

