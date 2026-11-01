











#ifndef SCATMECH_INSTRUMENT_H
#define SCATMECH_INSTRUMENT_H

#include "brdf.h"


namespace SCATMECH {



    class Instrument_BRDF_Model : public BRDF_Model
    {
        public:
            DECLARE_MODEL();
    };

    void Register(const Instrument_BRDF_Model* x);


}


#endif
