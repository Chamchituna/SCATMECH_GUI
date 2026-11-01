











#ifndef SCATMECH_LOCAL_H
#define SCATMECH_LOCAL_H

#include <complex>
#include "brdf.h"


namespace SCATMECH {


    class Local_BRDF_Model : public BRDF_Model {
        public:
            MuellerMatrix MuellerDSC(double thetai,
                                     double thetas,
                                     double phis,
                                     double rotation,
                                     Coordinate_System cs = psps);

            JonesMatrix JonesDSC(double thetai,
                                 double thetas,
                                 double phis,
                                 double rotation,
                                 Coordinate_System cs = psps);



            JonesMatrix JonesDSC(const Vector& source,
                                 const Vector& viewer,
                                 const Vector& normal,
                                 const Vector& xaxis,
                                 Coordinate_System cs=plane);



            MuellerMatrix MuellerDSC(const Vector& source,
                                     const Vector& viewer,
                                     const Vector& normal,
                                     const Vector& xaxis,
                                     Coordinate_System cs=plane);

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,density);

        protected:

            virtual void setup();

            virtual JonesMatrix jonesDSC() {
                error("Attempt to call jonesDSC for a depolarizing model.");
                return JonesMatrix();
            }

            virtual MuellerMatrix muellerDSC() {
                return this->jonesDSC();
            }

        private:

            JonesMatrix jones() {
                return this->jonesDSC()*sqrt(density/(fabs(cos(thetai)*cos(thetas))));
            }
            MuellerMatrix mueller() {
                return this->muellerDSC()*(density/(fabs(cos(thetai)*cos(thetas))));
            }
    };

    typedef Model_Ptr<Local_BRDF_Model> Local_BRDF_Model_Ptr;

    void Register(const Local_BRDF_Model* x);


}


#endif
