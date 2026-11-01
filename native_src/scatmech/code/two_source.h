












#ifndef TWO_SOURCE_H
#define TWO_SOURCE_H

#include "brdf.h"

namespace SCATMECH {

    class Two_Source_BRDF_Model : public BRDF_Model
    {
        public:
            MuellerMatrix mueller();

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,factor1);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source1);
            DECLARE_PARAMETER(double,factor2);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source2);
            DECLARE_PARAMETER(double,correlation);
    };

    class Three_Source_BRDF_Model : public BRDF_Model
    {
        public:
            MuellerMatrix mueller();

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,factor1);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source1);
            DECLARE_PARAMETER(double,factor2);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source2);
            DECLARE_PARAMETER(double,factor3);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source3);
    };

    class Four_Source_BRDF_Model : public BRDF_Model
    {
        public:
            MuellerMatrix mueller();

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,factor1);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source1);
            DECLARE_PARAMETER(double,factor2);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source2);
            DECLARE_PARAMETER(double,factor3);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source3);
            DECLARE_PARAMETER(double,factor4);
            DECLARE_PARAMETER(BRDF_Model_Ptr,source4);
    };

}

#endif

