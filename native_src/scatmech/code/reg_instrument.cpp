











#include "instrument.h"
#include "rayinst.h"
#include "finiteaperture.h"
#include "focussedbeam.h"


namespace SCATMECH {


    void Register(const Instrument_BRDF_Model* x)
    {
        static bool regd=false;
        if (!regd) {
            regd=true;

            Register_Model(Instrument_BRDF_Model);
            Register_Model(Rayleigh_Instrument_BRDF_Model);
            Register_Model(Finite_Aperture_Instrument_BRDF_Model);
            Register_Model(Focussed_Beam_Instrument_BRDF_Model);
        }
    }


}

