











#ifndef SCATMECH_OASPHERE_H
#define SCATMECH_OASPHERE_H

#include "sphrscat.h"
#include "dielfunc.h"
#include <vector>

namespace SCATMECH {









    class Optically_Active_Sphere_Scatterer: public Free_Space_Scatterer {
        public:
            Optically_Active_Sphere_Scatterer();

            virtual JonesMatrix jones(const Vector& kin,const Vector& kout);

            double CscaL();
            double CextL();
            double CbackL();
            double CabsL()
            {
                return CextL()-CscaL();
            }
            double CscaR();
            double CextR();
            double CbackR();
            double CabsR()
            {
                return CextR()-CscaR();
            }


            double QscaL() {
                return CscaL()/(pi*sqr(radius));
            }

            double QextL() {
                return CextL()/(pi*sqr(radius));
            }

            double QbackL() {
                return CbackL()/(pi*sqr(radius));
            }

            double QabsL() {
                return CabsL()/(pi*sqr(radius));
            }


            double QscaR() {
                return CscaR()/(pi*sqr(radius));
            }

            double QextR() {
                return CextL()/(pi*sqr(radius));
            }

            double QbackR() {
                return CbackR()/(pi*sqr(radius));
            }

            double QabsR() {
                return CabsR()/(pi*sqr(radius));
            }

        public:

            DECLARE_MODEL();

            DECLARE_PARAMETER(double,radius);

            DECLARE_PARAMETER(dielectric_function,materialleft);

            DECLARE_PARAMETER(dielectric_function,materialright);

        protected:

            std::vector<COMPLEX> a,b,c;
            std::vector<double> E,F;
            int NSTOP;
            int NMX;
            double k;
            double x;
            void setup();
    };

}


#endif
