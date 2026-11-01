











#ifndef SCATMECH_FLAKE_H
#define SCATMECH_FLAKE_H

#include "facet.h"


namespace SCATMECH {

    class Subsurface_Facet_BRDF_Model:
        public Facet_BRDF_Model
    {
        public:


            virtual double local_angle(double thetai,double thetas,double phis);
            virtual double local_slope(double thetai,double thetas,double phis);

            DECLARE_MODEL();


            DECLARE_PARAMETER(dielectric_function,overcoat);


            DECLARE_PARAMETER(StackModel_Ptr,overcoat_films);

        protected:


            void setup();


            JonesMatrix jones();
    };


}



#endif
