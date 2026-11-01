











#ifndef SCATMECH_NSPHERE_H
#define SCATMECH_NSPHERE_H

#include "sphrscat.h"
#include "dielfunc.h"
#include "filmtran.h"
#include <vector>

namespace SCATMECH {





    class MultilayerCoatedMieScatterer: public SphericalScatterer {
        public:
            MultilayerCoatedMieScatterer();


            COMPLEX s1(double angle);
            COMPLEX s2(double angle);

            double Csca();
            double Cext();
            double Cback();

        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(StackModel_Ptr,stack);

        private:

            std::vector<COMPLEX> A,B;
            std::vector<double> E,F;
            int NSTOP;
            int NMX;

            void setup();
    };
}


#endif
