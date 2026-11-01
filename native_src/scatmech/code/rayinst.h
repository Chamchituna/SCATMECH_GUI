











#ifndef SCATMECH_INSTSIG_H
#define SCATMECH_INSTSIG_H

#include "scatmech.h"
#include "instrument.h"
#include "filmtran.h"


namespace SCATMECH {






    class Rayleigh_Instrument_BRDF_Model : public Instrument_BRDF_Model
    {
        public:

            DECLARE_MODEL();


            DECLARE_PARAMETER(double,field_of_view);


            DECLARE_PARAMETER(dielectric_function,air);


            DECLARE_PARAMETER(double,number_density);


            DECLARE_PARAMETER(StackModel_Ptr,stack);

        protected:

            MuellerMatrix mueller();

        private:
    };


}



#endif
