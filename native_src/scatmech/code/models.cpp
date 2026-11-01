











#include "brdf.h"
#include "sphrscat.h"
#include "facet.h"
#include "rough.h"
#include "torrspar.h"
#include "rcw.h"
#include "crossrcw.h"
#include "reflectance.h"
#include "filmtran.h"

using namespace std;

namespace SCATMECH {

    void Register(const Model* x)
    {
        static int Models_Registered = 0;
        if (!Models_Registered) {
            Models_Registered=1;

            Register((BRDF_Model*)NULL);
            Register((Free_Space_Scatterer*)NULL);
            Register((Slope_Distribution_Function*)NULL);
            Register((PSD_Function*)NULL);
            Register((Shadow_Function*)NULL);
            Register((Reflectance*)NULL);
            Register((Grating*)NULL);
            Register((CrossGrating*)NULL);
            Register_Model(RCW_Model);
            Register_Model(CrossRCW_Model);
			Register((StackModel*)NULL);
        }
    }

}

