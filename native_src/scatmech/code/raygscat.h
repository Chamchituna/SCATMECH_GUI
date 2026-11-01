











#ifndef SCATMECH_RAYGSCAT_H
#define SCATMECH_RAYGSCAT_H

#include "scatmech.h"
#include "sphrscat.h"


namespace SCATMECH {



    class RayleighGansSphereScatterer: public SphericalScatterer {
        public:

            COMPLEX s1(double angle);
            COMPLEX s2(double angle);


            double Csca();
            double Cext();
            double Cback();

        public:

            DECLARE_MODEL();

        private:

            virtual void setup();
    };


}


#endif
