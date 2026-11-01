












#ifndef SCATMECH_AXIPART_H
#define SCATMECH_AXIPART_H

#include <stdio.h>
#include <vector>
#include "local.h"
#include "bobvlieg.h"
#include "tmatrix.h"


namespace SCATMECH {


    class Axisymmetric_Particle_BRDF_Model: public Local_BRDF_Model
    {
        public:

            COMPLEX Epp(double thetai,double thetas,double phis);
            COMPLEX Eps(double thetai,double thetas,double phis);
            COMPLEX Esp(double thetai,double thetas,double phis);
            COMPLEX Ess(double thetai,double thetas,double phis);

            DECLARE_MODEL();
            DECLARE_PARAMETER(Axisymmetric_Shape_Ptr,Shape);
            DECLARE_PARAMETER(dielectric_function,particle);
            DECLARE_PARAMETER(StackModel_Ptr,stack);
            DECLARE_PARAMETER(double,delta);
            DECLARE_PARAMETER(int,lmax);
            DECLARE_PARAMETER(int,mmax);
            DECLARE_PARAMETER(int,order);
            DECLARE_PARAMETER(int,Norm_Inc_Approx);
            DECLARE_PARAMETER(int,improve);

        public:
            Axisymmetric_Particle_BRDF_Model();

            MuellerMatrix Specular(double theta);




        protected:
            void setup();

            JonesMatrix jonesDSC();

        private:

            dielectric_function n3;

            COMPLEX N0;
            COMPLEX N1;
            COMPLEX N2;
            COMPLEX N3;


            double  old_thetai;
            double  old_thetas;
            double  old_phis;
            int     LMAX;
            int     MMAX;
            int     BH_LMAX;
            double  qq;
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
            ScatterTMatrix ScatMatrixInverse;
            std::vector<COMPLEX> Wp,Ws;
            std::vector<COMPLEX> Zp,Zs;
            std::vector<COMPLEX> Vp,Vs;
            std::vector<COMPLEX> eIP;




            COMPLEX VIp(int l,int m,int f,double thetai) const;
            COMPLEX VIRp(int l,int m,int f,double thetai) const;
            COMPLEX VITp(int l,int m,int f,double thetai) const;
            COMPLEX VIs(int l,int m,int f,double thetai) const;
            COMPLEX VIRs(int l,int m,int f,double thetai) const;
            COMPLEX VITs(int l,int m,int f,double thetai) const;

            void    Bmatrix(ScatterTMatrix& T);
            void    Amatrix(ScatterTMatrix& A);
            void    invert_block_diagonal(ScatterTMatrix& AA);
            void    set_geometry(double thetai,double thetas,double phis);

            COMPLEX E(std::vector<COMPLEX>& W,std::vector<COMPLEX>& Z);
            void    calculate_W(double thetai);
            void    calculate_Z(double thetas);
            void    calculate_eIP(double phis);


            void    iterative_improvement(ScatterTMatrix& Ainv, ScatterTMatrix& A,
                                          std::vector<COMPLEX>& b,std::vector<COMPLEX>& x);

    };


}



#endif
