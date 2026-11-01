











#ifndef SCATMECH_RAYSCAT_H
#define SCATMECH_RAYSCAT_H

#include "scatmech.h"
#include "sphrscat.h"


namespace SCATMECH {






    class RayleighScatterer: public SphericalScatterer {
        public:


            COMPLEX s1(double angle);
            COMPLEX s2(double angle);


            double Csca();
            double Cext();
            double Cback();

        public:
            DECLARE_MODEL();

        private:

            COMPLEX polarizability;

            virtual void setup();
    };


}


#endif
