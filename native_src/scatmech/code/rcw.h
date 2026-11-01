











#ifndef SCATMECH_RCW_H
#define SCATMECH_RCW_H

#include "grating.h"
#include "brdf.h"
#include <vector>
#include <cfloat>

namespace SCATMECH {

























    class RCW_Model : public Model {
        public:

            JonesMatrix GetAmplitude(int i) {
                SETUP();
                if (i < -order || i > order) return JonesZero();
                return r[i+order];
            }

            MuellerMatrix GetIntensity(int i) {
                SETUP();
                if (i < -order || i > order) return MuellerZero();
                return R[i+order];
            }

            int GetMinimumPropagatingOrder() {
                SETUP();
                return min_order;
            }

            int GetMaximumPropagatingOrder() {
                SETUP();
                return max_order;
            }

            Vector GetDirection(int i) {
                SETUP();
                if (i < -order || i > order) return Vector(0.,0.,0.);
                return V[i+order];
            }

            CVector GetPropagationVector(int i) {
                SETUP();
                if (i < -order || i > order) return CVector(0.,0.,0.);
                return kvector[i+order];
            }

            CVector GetEField(const JonesVector& inpol, const Vector& r, bool incident=false);
            CVector GetHField(const JonesVector& inpol, const Vector& r, bool incident=false);
            CVector GetBField(const JonesVector& inpol, const Vector& r, bool incident=false);
            CVector GetDField(const JonesVector& inpol, const Vector& r, bool incident=false);

            DECLARE_MODEL();
            DECLARE_PARAMETER(int,order);
            DECLARE_PARAMETER(int,type);
            DECLARE_PARAMETER(double,lambda);
            DECLARE_PARAMETER(Grating_Ptr,grating);
            DECLARE_PARAMETER(double,thetai);
            DECLARE_PARAMETER(double,rotation);

        protected:
            void setup();

        private:
            std::vector<JonesMatrix> r;
            std::vector<MuellerMatrix> R;
            std::vector<Vector> V;
            std::vector<CVector> kvector;

            double k0,inckx,incky,inckz;
            int n;
            int nmat;
            double period;
            int min_order;
            int max_order;
            double totalthickness;

        private:
            void InPlaneReflection(bool backward);
            void InPlaneTransmission(bool backward);
            void ConicalReflection(bool backward);
            void ConicalTransmission(bool backward);
            void ConicalAnisoReflection(bool backward);
            void ConicalAnisoTransmission(bool backward);

        private:
            double nI,eI;
            COMPLEX nII,eII;

            std::vector<double> kxi,Kx;
            std::vector<COMPLEX> kIzi,kIIzi,YI,YII,ZI,ZII;
            double ky;
    };

    class RCW_BRDF_Model: public BRDF_Model {
        public:
            RCW_BRDF_Model() {
                set_recalc_on_direction_change(0xF0,0x00);
                RCW.init();
            }

            MuellerMatrix mueller();

        protected:

            void setup();

        protected:
            RCW_Model RCW;

            double cosalpha;
            double Omega;

            int low,high;

            std::vector<MuellerMatrix> R;
            std::vector<Vector> V;

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,alpha);
            DECLARE_PARAMETER(int,order);
            DECLARE_PARAMETER(Grating_Ptr,grating);
    };

}

#endif
