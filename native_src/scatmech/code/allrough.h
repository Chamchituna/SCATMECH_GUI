











#ifndef SCATMECH_GRADEDROUGH_H
#define SCATMECH_GRADEDROUGH_H

#include "roughnes.h"


namespace SCATMECH {


    class Correlated_Roughness_Stack_BRDF_Model : public Roughness_BRDF_Model
    {
        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(StackModel_Ptr,stack);


        protected:

            void setup();
            JonesMatrix jones();

        private:

            Roughness_Stack_BRDF_Model model;
    };

    class Uncorrelated_Roughness_Stack_BRDF_Model : public Roughness_BRDF_Model
    {
        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(StackModel_Ptr,stack);

        protected:

            void setup();
            MuellerMatrix mueller();

        private:

            Roughness_Stack_BRDF_Model model;
    };

    class Growth_Roughness_Stack_BRDF_Model : public Roughness_BRDF_Model
    {
        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(StackModel_Ptr,stack);
            DECLARE_PARAMETER(double,relaxation);
            DECLARE_PARAMETER(double,exponent);
            DECLARE_PARAMETER(PSD_Function_Ptr,intrinsic);

        protected:

            void setup();
            MuellerMatrix mueller();

        private:

            Roughness_Stack_BRDF_Model model;
    };

}


#endif
