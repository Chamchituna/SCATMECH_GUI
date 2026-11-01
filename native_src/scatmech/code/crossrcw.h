











#ifndef CROSSRCW_H
#define CROSSRCW_H

#include "crossgrating.h"
#include "mueller.h"
#include "vector3d.h"
#include "brdf.h"
#include <limits>

namespace SCATMECH {









    class CrossRCW_Model : public Model {
        public:
            CrossRCW_Model() : old_type(-1) {}

            JonesMatrix GetAmplitude(int i,int j);
            MuellerMatrix GetIntensity(int i,int j);
			StokesVector GetAbsorption();
            Vector GetDirection(int i,int j);
            CVector GetPropagationVector(int i,int j);

            CVector GetEField(const JonesVector& inpol, const Vector& r, bool incident=false);
            CVector GetHField(const JonesVector& inpol, const Vector& r, bool incident=false);
            CVector GetBField(const JonesVector& inpol, const Vector& r, bool incident=false);
            CVector GetDField(const JonesVector& inpol, const Vector& r, bool incident=false);

        protected:
            void setup();

            virtual void set_parameter_base(
                const STRING& parameter,
                const STRING& value
            );

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,thetai);
            DECLARE_PARAMETER(double,rotation);
            DECLARE_PARAMETER(double,lambda);
            DECLARE_PARAMETER(int,type);
            DECLARE_PARAMETER(int,order1);
            DECLARE_PARAMETER(int,order2);
            DECLARE_PARAMETER(CrossGrating_Ptr,grating);
        private:

            std::string write_eigenvalues;
            int old_type;
            double totalthickness;
            double nI,eI;
            COMPLEX nII,eII;
            double k0;


            FARRAY<CVector> Vr;
            FARRAY<CVector> Vt;
            FARRAY<MuellerMatrix> R;
            FARRAY<JonesMatrix> r;
            FARRAY<MuellerMatrix> T;
            FARRAY<JonesMatrix> t;

    };

    typedef Model_Ptr<CrossRCW_Model> CrossRCW_Model_Ptr;

    class CrossRCW_BRDF_Model: public BRDF_Model {
        public:
            CrossRCW_BRDF_Model() {
                set_recalc_on_direction_change(0xF0,0x00);
            }

            MuellerMatrix mueller();

        protected:

            void setup();

        protected:
            double cosalpha;
            double Omega;

            CrossRCW_Model RCW;

            int M1,M2;

            FARRAY<MuellerMatrix> R;
            FARRAY<Vector> V;

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,alpha);
            DECLARE_PARAMETER(int,order1);
            DECLARE_PARAMETER(int,order2);
            DECLARE_PARAMETER(CrossGrating_Ptr,grating);
    };
}

#endif
