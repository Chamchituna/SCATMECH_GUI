











#ifndef SCATMECH_SUBBOBVLIEG_H
#define SCATMECH_SUBBOBVLIEG_H

#include "bobvlieg.h"
#include "axipart.h"

namespace SCATMECH {

    class Subsurface_Bobbert_Vlieger_BRDF_Model : public Local_BRDF_Model
    {
        public:
            DECLARE_MODEL();
            DECLARE_PARAMETER(dielectric_function,sphere);
            DECLARE_PARAMETER(double,radius);
            DECLARE_PARAMETER(StackModel_Ptr,spherecoat);
            DECLARE_PARAMETER(StackModel_Ptr,stack);
            DECLARE_PARAMETER(double,delta);
            DECLARE_PARAMETER(int,lmax);
            DECLARE_PARAMETER(int,order);
            DECLARE_PARAMETER(int,Norm_Inc_Approx);
            DECLARE_PARAMETER(int,improve);

        public:

            MuellerMatrix Specular(double theta) {
                return model.Specular(theta);
            }

        protected:
            void setup();
            JonesMatrix jonesDSC();

        private:
            Bobbert_Vlieger_BRDF_Model model;
    };

    class Subsurface_Axisymmetric_Particle_BRDF_Model : public Local_BRDF_Model
    {
        public :
            DECLARE_MODEL();
            DECLARE_PARAMETER(Axisymmetric_Shape_Ptr,Shape);
            DECLARE_PARAMETER(dielectric_function,particle);
            DECLARE_PARAMETER(StackModel_Ptr,stack);
            DECLARE_PARAMETER(double,delta);
            DECLARE_PARAMETER(int,lmax);
            DECLARE_PARAMETER(int,mmax);
            DECLARE_PARAMETER(int,order);
            DECLARE_PARAMETER(int,Norm_Inc_Approx);
            DECLARE_PARAMETER(int,improve);

        public:

            MuellerMatrix Specular(double theta) {
                return model.Specular(theta);
            }

        protected:
            void setup();
            JonesMatrix jonesDSC();

        private:
            Axisymmetric_Particle_BRDF_Model model;
    };
}

#endif
