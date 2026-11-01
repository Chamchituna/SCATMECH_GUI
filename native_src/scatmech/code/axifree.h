












#ifndef SCATMECH_AXIFREE_H
#define SCATMECH_AXIFREE_H

#include <stdio.h>
#include <vector>
#include "tmatrix.h"
#include "scatmatrix.h"


namespace SCATMECH {

    class TMatrix_Axisymmetric_Scatterer : public Free_Space_Scatterer
    {
        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(Axisymmetric_Shape_Ptr,Shape);
            DECLARE_PARAMETER(dielectric_function,particle);
            DECLARE_PARAMETER(int,lmax);
            DECLARE_PARAMETER(int,mmax);

        public:
            TMatrix_Axisymmetric_Scatterer();

            virtual JonesMatrix jones(const Vector& kin,const Vector& kout);

        protected:
            void setup();

        private:

            COMPLEX N2;
            double N0;

            double  a;


            double  old_thetai;
            double  old_thetas;
            double  old_phis;

            int     LMAX;
            int     MMAX;
            int     BH_LMAX;
            double  k;
            int     sqrsize;
            int     old_LMAX;


            int index(int l,int m, int f) const {
                return 2*(l*l+l-1+m)+f;
            }
            int index(int l,int m) const {
                return l*l+l-1+m;
            }




            ScatterTMatrix ScatMatrix;
            std::vector<COMPLEX> Wp,Ws;
            std::vector<COMPLEX> Zp,Zs;
            std::vector<COMPLEX> Vp,Vs;
            std::vector<COMPLEX> eIP;

            std::vector<double> lvector;


            COMPLEX VIp(int l,int m,int f,double thetai) const;
            COMPLEX VIs(int l,int m,int f,double thetai) const;

            void    Bmatrix(ScatterTMatrix& T);

            void    set_geometry(double thetai,double thetas,double phis);

            COMPLEX E(std::vector<COMPLEX>& W,std::vector<COMPLEX>& Z);
            void    calculate_W(double thetai);
            void    calculate_Z(double thetas);
            void    calculate_eIP(double phis);
    };


}



#endif
