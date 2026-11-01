











#include <float.h>
#include <assert.h>
#include <vector>
#include <cstdlib>
#include "scatmech.h"
#include "axipart.h"
#include "bobvlieg.h"

using namespace std;


namespace SCATMECH {

    using namespace BobVlieg_Supp;






    COMPLEX
    Axisymmetric_Particle_BRDF_Model::
    VIp(int l,int m,int f,double _thetai) const
    {
        COMPLEX result;
        COMPLEX costheta2=cos(_thetai/2.);
        COMPLEX sintheta2=sin(_thetai/2.);
        if (f==efield) {

            result =  1./k*ipow(l-1)*lvector[l]*mpow(m-1)*
                      dminus(l,m,costheta2,sintheta2);
        } else {

            result =  -1./k*ipow(l)*lvector[l]*mpow(m-1)*
                      dplus(l,m,costheta2,sintheta2);
        }
        return result;
    }





    COMPLEX
    Axisymmetric_Particle_BRDF_Model::
    VIRp(int l,int m,int f,double _thetai) const
    {
        COMPLEX result;
        COMPLEX cost = cos(_thetai);
        COMPLEX phase = exp(2.*cI*qq*cost);
        COMPLEX rp = stack->rp12(_thetai,lambda,vacuum,substrate);
        COMPLEX costheta2=cos(_thetai/2.);
        COMPLEX sintheta2=sin(_thetai/2.);

        if (f==efield) {

            result = 1./k*phase*rp*
                     ipow(l-1)*lvector[l]*mpow(l)*
                     dminus(l,m,costheta2,sintheta2);
        } else {

            result =  -1./k*phase*rp*
                      ipow(l)*lvector[l]*mpow(l-1)*
                      dplus(l,m,costheta2,sintheta2);
        }
        return result;
    }





    COMPLEX
    Axisymmetric_Particle_BRDF_Model::
    VITp(int l,int m,int f,double _thetai) const
    {
        COMPLEX sint = sin(_thetai)*substrate.n(lambda);
        COMPLEX cost = sqrt(1.-sqr(sint));
        if (imag(cost)<0) cost = -cost;
        COMPLEX __thetai = arcsine(sint);
        if (imag(__thetai)>0) __thetai = conj(__thetai);
        COMPLEX sintheta2=sqrt(1.-cost)/sqrt(2.);
        if (real(sintheta2)<0) sintheta2 = -sintheta2;
        COMPLEX costheta2=sint/2./sintheta2;
        COMPLEX _tp = stack->tp12(__thetai,lambda,vacuum,substrate);
        COMPLEX phase = exp(cI*qq*(cost-substrate.n(lambda)*cos(_thetai)));

        COMPLEX result;
        if (f==efield) {

            result =  _tp/k*
                      ipow(l-1)*lvector[l]*mpow(l)*
                      dminus(l,m,costheta2,sintheta2);
        } else {

            result =  -_tp/k*
                      ipow(l)*lvector[l]*mpow(l-1)*
                      dplus(l,m,costheta2,sintheta2);
        }
        return result*phase;
    }





    COMPLEX
    Axisymmetric_Particle_BRDF_Model::
    VIs(int l,int m,int f,double _thetai) const
    {
        COMPLEX result;
        COMPLEX costheta2=cos(_thetai/2.);
        COMPLEX sintheta2=sin(_thetai/2.);

        if (f==efield) {

            result =  -cI/k*
                      ipow(l-1)*lvector[l]*mpow(m-1)*
                      dplus(l,m,costheta2,sintheta2);
        } else {

            result = cI/k*
                     ipow(l)*lvector[l]*mpow(m-1)*
                     dminus(l,m,costheta2,sintheta2);
        }
        return result;
    }





    COMPLEX
    Axisymmetric_Particle_BRDF_Model::
    VIRs(int l,int m,int f,double _thetai) const
    {
        COMPLEX result;
        COMPLEX cost = cos(_thetai);
        COMPLEX phase = exp(2.*cI*qq*cost);
        COMPLEX rs = stack->rs12(_thetai,lambda,vacuum,substrate);
        COMPLEX costheta2=cos(_thetai/2.);
        COMPLEX sintheta2=sin(_thetai/2.);

        if (f==efield) {

            result =  -cI/k*phase*rs*
                      ipow(l-1)*lvector[l]*
                      mpow(l-1)*dplus(l,m,costheta2,sintheta2);
        } else {

            result =  cI/k*phase*rs*
                      ipow(l)*lvector[l]*mpow(l)*
                      dminus(l,m,costheta2,sintheta2);
        }
        return result;
    }





    COMPLEX
    Axisymmetric_Particle_BRDF_Model::
    VITs(int l,int m,int f,double _thetai) const
    {
        COMPLEX sint = sin(_thetai)*substrate.n(lambda);
        COMPLEX cost = sqrt(1.-sqr(sint));
        if (imag(cost)<0) cost = -cost;
        COMPLEX __thetai = arcsine(sint);
        if (imag(__thetai)>0) __thetai = conj(__thetai);
        COMPLEX sintheta2=sqrt(1.-cost)/sqrt(2.);
        if (real(sintheta2)<0) sintheta2 = -sintheta2;
        COMPLEX costheta2=sint/2./sintheta2;
        COMPLEX _ts = stack->ts12(__thetai,lambda,vacuum,substrate);
        COMPLEX phase = exp(cI*qq*(cost-substrate.n(lambda)*cos(_thetai)));

        COMPLEX result;
        if (f==efield) {

            result =  -cI/k*_ts*
                      ipow(l-1)*lvector[l]*
                      mpow(l-1)*dplus(l,m,costheta2,sintheta2);
        } else {

            result = cI/k*_ts*
                     ipow(l)*lvector[l]*mpow(l)*
                     dminus(l,m,costheta2,sintheta2);
        }
        return result*phase;
    }




    void
    Axisymmetric_Particle_BRDF_Model::
    iterative_improvement(ScatterTMatrix& Ainv, ScatterTMatrix& A, vector<COMPLEX>& b, vector<COMPLEX>& x)
    {
        vector<COMPLEX> temp1(sqrsize);
        int m,l,f,l_,f_;

        for (m=-MMAX; m<=MMAX; ++m) {
            COMPLEX temp2;
            int beginl = (m==0) ? 1 : abs(m);
            for (l_=beginl; l_<=LMAX; ++l_) {
                for (f_=0; f_<=1; ++f_) {
                    int i_ = index(l_,m,f_);
                    temp1[i_]=0;
                    for (l=beginl; l<=LMAX; ++l) {
                        for (f=0; f<=1; ++f) {
                            int i = index(l,m,f);
                            int mm = LMAX+m;
                            int ll = 2*(LMAX-l)+f;
                            int ll_ = 2*(LMAX-l_)+f_;

                            temp1[i_]+=A[mm][ll_][ll]*x[i];
                        }
                    }
                    temp1[i_]-=b[i_];
                }
            }
            for (l_=beginl; l_<=LMAX; ++l_) {
                for (f_=0; f_<=1; ++f_) {
                    int i_ = index(l_,m,f_);
                    temp2=0;
                    for (l=beginl; l<=LMAX; ++l) {
                        for (f=0; f<=1; ++f) {
                            int i = index(l,m,f);
                            int mm = LMAX+m;
                            int ll = 2*(LMAX-l)+f;
                            int ll_ = 2*(LMAX-l_)+f_;

                            temp2+=Ainv[mm][ll_][ll]*temp1[i];
                        }
                    }
                    x[i_]-=temp2;
                }
            }
        }
    }





    void
    Axisymmetric_Particle_BRDF_Model::
    calculate_W(double _thetai)
    {



        if (is_down_to_up()||is_down_to_down()) {
            for (int l=1; l<=LMAX; ++l) {
                for (int m=-l; m<=l; ++m) {
                    for (int f=0; f<=1; ++f) {
                        int i=index(l,m,f);
                        Vp[i]=VIp(l,m,f,_thetai)+VIRp(l,m,f,_thetai);
                        Vs[i]=VIs(l,m,f,_thetai)+VIRs(l,m,f,_thetai);
                    }
                }
            }
        } else {
            COMPLEX sint = sin(_thetai)*substrate.n(lambda);
            COMPLEX cost = sqrt(1.-sqr(sint));



            double factor = abs(COMPLEX(cos(_thetai))/cost*COMPLEX(sqrt(substrate.n(lambda))));

            for (int l=1; l<=LMAX; ++l) {
                for (int m=-l; m<=l; ++m) {
                    for (int f=0; f<=1; ++f) {
                        int i=index(l,m,f);
                        Vp[i]=VITp(l,m,f,_thetai)*factor;
                        Vs[i]=VITs(l,m,f,_thetai)*factor;
                    }
                }
            }
        }



        for (int m=-MMAX; m<=MMAX; ++m) {
            int beginl = (m==0) ? 1 : abs(m);
            for (int l_=beginl; l_<=LMAX; ++l_) {
                for (int f_=0; f_<=1; ++f_) {
                    int i_ = index(l_,m,f_);
                    Wp[i_]=0;
                    Ws[i_]=0;
                    for (int l=beginl; l<=LMAX; ++l) {
                        for (int f=0; f<=1; ++f) {
                            int i = index(l,m,f);
                            int mm = LMAX+m;
                            int ll = 2*(LMAX-l)+f;
                            int ll_ = 2*(LMAX-l_)+f_;

                            Wp[i_]+=ScatMatrix[mm][ll_][ll]*Vp[i];
                            Ws[i_]+=ScatMatrix[mm][ll_][ll]*Vs[i];
                        }
                    }
                }
            }
        }



        if (order<0) {
            for (int i=0; i<improve; ++i) {
                iterative_improvement(ScatMatrix,ScatMatrixInverse,Vp,Wp);
                iterative_improvement(ScatMatrix,ScatMatrixInverse,Vs,Ws);
            }
        }
    }






    void
    Axisymmetric_Particle_BRDF_Model::
    calculate_Z(double _thetas)
    {
        if (is_down_to_up()||is_up_to_up()) {
            double d_cost = cos(_thetas);
            COMPLEX cost= d_cost;
            double sint = real(sqrt(1.-sqr(cost)));

            COMPLEX _cost=cos(pi-_thetas);
            COMPLEX rp = stack->rp12(pi-_thetas,lambda,vacuum,substrate);
            COMPLEX rs = stack->rs12(pi-_thetas,lambda,vacuum,substrate);
            COMPLEX phase = exp(2.*cI*qq*_cost);
            COMPLEX rpphase = rp*phase;
            COMPLEX rsphase = rs*phase;

            for (int i=0,l=1; l<=LMAX; ++l) {
                if (l<=BH_LMAX) {
                    int _mmax = (MMAX>l) ? l : MMAX;
                    for (int m=-_mmax; m<=_mmax; ++m) {
                        COMPLEX temp1 = Ptilde(l,m,cost)*(double)m/sint;
                        COMPLEX temp2 = P_tilde(l,m,cost);
                        for (int f=0; f<=1; ++f,++i) {
                            if (f==efield) {
                                Zp[i] = ipow(l)*mpow(l)*temp2;
                                Zs[i] = -ipow(l+1)*mpow(l+1)*temp1;
                                Zp[i] = Zp[i]*(1.+mpow(l-m+1)*rpphase);
                                Zs[i] = Zs[i]*(1.+mpow(l-m)*rsphase);
                            } else {
                                Zp[i] = -ipow(l+1)*mpow(l+1)*temp1;
                                Zs[i] = -ipow(l)*mpow(l)*temp2;
                                Zp[i] = Zp[i]*(1.+mpow(l-m)*rpphase);
                                Zs[i] = Zs[i]*(1.+mpow(l-m+1)*rsphase);
                            }
                        }
                    }
                } else {
                    int _mmax = (MMAX>l) ? l : MMAX;
                    for (int m=-_mmax; m<=_mmax; ++m) {
                        for (int f=0; f<=1; ++f,++i) {
                            Zp[i] = 0.;
                            Zs[i] = 0.;
                            Zp[i] = 0.;
                            Zs[i] = 0.;
                        }
                    }
                }
            }
        } else {
            double index = substrate.n(lambda);
            COMPLEX sint = (COMPLEX)(sin(pi-_thetas)*index);
            COMPLEX __thetas = arcsine(sint);
            if (imag(__thetas)>0) __thetas = conj(__thetas);
            COMPLEX cost= sqrt(1.-sqr(sint));
            if (imag(cost)<0) cost = -cost;
            COMPLEX tp = stack->tp12(__thetas,lambda,vacuum,substrate);
            COMPLEX ts = stack->ts12(__thetas,lambda,vacuum,substrate);
            COMPLEX phase = exp(cI*qq*(cost-substrate.n(lambda)*cos(pi-_thetas)));


            double factor = abs(COMPLEX(cos(pi-_thetas))/cost*COMPLEX(sqrt(cube(index))));
            COMPLEX tpphase = tp*phase*factor;
            COMPLEX tsphase = ts*phase*factor;

            for (int i=0,l=1; l<=LMAX; ++l) {
                if (l<=BH_LMAX) {
                    for (int m=-l; m<=l; ++m) {
                        COMPLEX temp1 = Ptilde(l,m,cost)*(double)m/sint;
                        COMPLEX temp2 = P_tilde(l,m,cost);
                        for (int f=0; f<=1; ++f,++i) {
                            if (f==efield) {
                                Zp[i] = ipow(l)*mpow(l)*temp2*tpphase;
                                Zs[i] = -ipow(l+1)*mpow(l+1)*temp1*tsphase;
                            } else {
                                Zp[i] = -ipow(l+1)*mpow(l+1)*temp1*tpphase;
                                Zs[i] = -ipow(l)*mpow(l)*temp2*tsphase;
                            }
                        }
                    }
                } else {
                    for (int m=-l; m<=l; ++m) {
                        for (int f=0; f<=1; ++f,++i) {
                            Zp[i] = 0.;
                            Zs[i] = 0.;
                            Zp[i] = 0.;
                            Zs[i] = 0.;
                        }
                    }
                }
            }
        }

    }








    void
    Axisymmetric_Particle_BRDF_Model::
    calculate_eIP(double _phis)
    {
        vector<COMPLEX> expmphi(2*LMAX+2);

        COMPLEX expphi=exp(cI*_phis);
        expmphi[LMAX]=1.;

        for (int m=1; m<=LMAX; ++m) {
            expmphi[LMAX+m]=expmphi[LMAX+m-1]*expphi;
            expmphi[LMAX-m]=expmphi[LMAX-m+1]/expphi;
        }

        for (int i=0,l=1; l<=BH_LMAX; ++l) {
            int _mmax = (MMAX>l) ? l : MMAX;
            for (int m=-_mmax; m<=_mmax; ++m) {
                for (int f=0; f<=1; ++f,++i) {
                    eIP[i] = expmphi[m+LMAX];
                }
            }
        }
    }





    COMPLEX
    Axisymmetric_Particle_BRDF_Model::
    E(vector<COMPLEX>& W,vector<COMPLEX>& Z)
    {
        COMPLEX E(0.);

        for (int i=0,l=1; l<=BH_LMAX; ++l) {
            int _mmax = (MMAX>l) ? l : MMAX;
            for (int m=-_mmax; m<=_mmax; ++m) {
                for (int f=0; f<=1; ++f,++i) {
                    E += W[i]*Z[i]*eIP[i];
                }
            }
        }

        return COMPLEX(0, -1) * E;
    }


}




