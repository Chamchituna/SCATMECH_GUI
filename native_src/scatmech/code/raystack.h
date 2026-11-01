











#ifndef SCATMECH_STACK_BR_H
#define SCATMECH_STACK_BR_H

#include "scatmech.h"
#include "filmtran.h"
#include "local.h"

namespace SCATMECH {




    class Rayleigh_Stack_BRDF_Model : public Local_BRDF_Model
    {
        public:

            virtual void setup();


            SCATMECH::JonesMatrix jonesDSC();

            DECLARE_MODEL();
            DECLARE_PARAMETER(StackModel_Ptr,stack);
            DECLARE_PARAMETER(double,depth);
            DECLARE_PARAMETER(double,radius);
            DECLARE_PARAMETER(dielectric_function,sphere);

        protected:



            double d;
            double lambda_eff;


            dielectric_function material1;
            dielectric_function material2;
            dielectric_function material3;
            dielectric_function materialS;



            double tau;


            dielectric_stack below;


            dielectric_stack above;


    };

}

#endif
