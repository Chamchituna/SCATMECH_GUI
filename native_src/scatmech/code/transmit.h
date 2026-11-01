














#ifndef SCATMECH_TRANSMIT_H

#include "brdf.h"
#include "filmtran.h"

namespace SCATMECH {

    class Transmit_BRDF_Model : public BRDF_Model
    {
        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(BRDF_Model_Ptr,model);
            DECLARE_PARAMETER(StackModel_Ptr,films);

        protected:
            virtual MuellerMatrix mueller();
    };

}

#endif

