











#ifndef SCATMECH_LAMBERT_H
#define SCATMECH_LAMBERT_H

#include "brdf.h"
#include "scattabl.h"
#include "reflectance.h"

namespace SCATMECH {



    class Lambertian_BRDF_Model : public BRDF_Model {
        public:
            DECLARE_MODEL();



            DECLARE_PARAMETER(Reflectance_Ptr,reflectance);
        protected:

            MuellerMatrix mueller();
    };

    typedef Model_Ptr<Lambertian_BRDF_Model> Lambertian_BRDF_Model_Ptr;

    void Register(const Lambertian_BRDF_Model* x);


}


#endif
