











#ifndef SCATMECH_COATED_H
#define SCATMECH_COATED_H

#include "sphrscat.h"
#include "dielfunc.h"
#include <vector>


namespace SCATMECH {







    class CoatedMieScatterer: public SphericalScatterer {
        public:
            CoatedMieScatterer();


            COMPLEX s1(double angle);
            COMPLEX s2(double angle);

            double Csca();
            double Cext();
            double Cback();

        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,thickness);
            DECLARE_PARAMETER(dielectric_function,coating);
            DECLARE_PARAMETER(int,nmax);

        private:

            std::vector<COMPLEX> A,B;
            std::vector<double> E,F;
            int NSTOP;
            int NMX;

            void setup();
    };


}


#endif
