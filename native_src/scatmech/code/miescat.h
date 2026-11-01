











#ifndef SCATMECH_MIESCAT_H
#define SCATMECH_MIESCAT_H

#include "sphrscat.h"
#include "dielfunc.h"
#include <vector>


namespace SCATMECH {







    class MieScatterer: public SphericalScatterer {
        public:

            MieScatterer();


            COMPLEX s1(double angle);
            COMPLEX s2(double angle);

            double Csca();
            double Cext();
            double Cback();

        public:

            DECLARE_MODEL();

        private:

            std::vector<COMPLEX> A,B;
            std::vector<double> E,F;
            int NSTOP;
            int NMX;

            virtual void setup();
    };


}


#endif
