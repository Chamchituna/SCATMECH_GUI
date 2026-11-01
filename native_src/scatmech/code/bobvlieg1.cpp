












#include <fstream>
#include <vector>
#include <sstream>
#include <cstdlib>
#include "scatmech.h"
#include "bobvlieg.h"
#include "dielfunc.h"
#include "askuser.h"

using namespace std;


namespace SCATMECH {

    using namespace BobVlieg_Supp;












    COMPLEX
    Bobbert_Vlieger_BRDF_Model::
    Bmatrix(int l_,int m_,int f_,int l,int m,int f) const
    {
        COMPLEX result;
        if (l_!=l || m_!=m || f_!=f) {
            result = 0.;
        } else if (d==0.) {
            COMPLEX psi_lntildeq = psi_(l,N2*q);
            COMPLEX psilntildeq = psi(l,N2*q);
            COMPLEX psi_lq = psi_(l,q);
            COMPLEX psilq = psi(l,q);
            COMPLEX zetalq = zeta(l,q);
            COMPLEX zeta_lq = zeta_(l,q);

            if (f==efield) {
                result =  -(N2 * psi_lq*psilntildeq -psilq*psi_lntildeq)/
                          (N2 * zeta_lq*psilntildeq-zetalq*psi_lntildeq);
            } else {
                result =  -(N2 * psilq*psi_lntildeq -psi_lq*psilntildeq)/
                          (N2 * zetalq*psi_lntildeq-zeta_lq*psilntildeq);
            }
        } else {







            int n = l;
            int rr = spherecoat->get_n()+1;
            if (f==efield) {
                COMPLEX T = 0.;
                double Rs = r0;
                COMPLEX ns = N2;
                for (int s=1; s<=rr-1; ++s) {
                    COMPLEX nsp1 = spherecoat->get_e()[s-1].index(lambda);
                    COMPLEX ms = nsp1/ns;
                    COMPLEX ys = 2*pi*ns/lambda*Rs;

                    T = -(ms*psi(n,ms*ys)*(psi_(n,ys)+T*chi_(n,ys))-psi_(n,ms*ys)*(psi(n,ys)+T*chi(n,ys)))/
                        (ms*chi(n,ms*ys)*(psi_(n,ys)+T*chi_(n,ys))-chi_(n,ms*ys)*(psi(n,ys)+T*chi(n,ys)));

                    ns = nsp1;
                    Rs += spherecoat->get_t()[s-1];
                }
                COMPLEX nr = spherecoat->get_e()[rr-2].index(lambda);
                COMPLEX mr = 1./nr;
                COMPLEX yr = 2*pi*nr/lambda*Rs;

                COMPLEX an = (mr*psi(n,mr*yr)*(psi_(n,yr)+T*chi_(n,yr))-psi_(n,mr*yr)*(psi(n,yr)+T*chi(n,yr)))/
                             (mr*zeta(n,mr*yr)*(psi_(n,yr)+T*chi_(n,yr))-zeta_(n,mr*yr)*(psi(n,yr)+T*chi(n,yr)));
                result = -an;
            } else if (f==hfield) {
                COMPLEX S = 0.;
                double Rs = r0;
                COMPLEX ns = N2;
                for (int s=1; s<=rr-1; ++s) {
                    COMPLEX nsp1 = spherecoat->get_e()[s-1].index(lambda);
                    COMPLEX ms = nsp1/ns;
                    COMPLEX ys = 2*pi*ns/lambda*Rs;

                    S = -(psi(n,ms*ys)*(psi_(n,ys)+S*chi_(n,ys))-ms*psi_(n,ms*ys)*(psi(n,ys)+S*chi(n,ys)))/
                        (chi(n,ms*ys)*(psi_(n,ys)+S*chi_(n,ys))-ms*chi_(n,ms*ys)*(psi(n,ys)+S*chi(n,ys)));
                    ns = nsp1;
                    Rs += spherecoat->get_t()[s-1];
                }
                COMPLEX nr = spherecoat->get_e()[rr-2].index(lambda);
                COMPLEX mr = 1./nr;
                COMPLEX yr = 2*pi*nr/lambda*Rs;


                COMPLEX bn = (psi(n,mr*yr)*(psi_(n,yr)+S*chi_(n,yr))-mr*psi_(n,mr*yr)*(psi(n,yr)+S*chi(n,yr)))/
                             (zeta(n,mr*yr)*(psi_(n,yr)+S*chi_(n,yr))-mr*zeta_(n,mr*yr)*(psi(n,yr)+S*chi(n,yr)));
                result = -bn;
            }
        }

        return result;
    }




    void
    Bobbert_Vlieger_BRDF_Model::
    Amatrix(ScatterTMatrix& _A)
    {


        int npts=((LMAX/10)+1)*10;
        if (npts>100) npts=100;
        int ipt = npts/10-1;



        vector<COMPLEX> U1(index(LMAX,LMAX)+1);
        vector<COMPLEX> U2(index(LMAX,LMAX)+1);
        vector<COMPLEX> U3(index(LMAX,LMAX)+1);
        vector<COMPLEX> U4(index(LMAX,LMAX)+1);
        vector<COMPLEX> V1(index(LMAX,LMAX)+1);
        vector<COMPLEX> V2(index(LMAX,LMAX)+1);
        for (int l=1; l<=LMAX; ++l) {
            for (int m=0; m<=l; ++m) {
                int i=index(l,m);

                double temp1=(2.*l-1.)*(2.*l+1.);
                double temp2=(2.*l+3.)*(2.*l+1.);

                U1[i] = (l-1.)/2.*ipow(-(l-1))*sqrt((l+m-1.)*(l+m)/temp1);
                U2[i] = (l-1.)/2.*ipow(-(l-1))*sqrt((l-m-1.)*(l-m)/temp1);
                U3[i] =-(l+2.)/2.*ipow(-(l+1))*sqrt((l-m+1.)*(l-m+2.)/temp2);
                U4[i] =-(l+2.)/2.*ipow(-(l+1))*sqrt((l+m+1.)*(l+m+2.)/temp2);

                V1[i] = 0.5*ipow(-(l-1))*sqrt((l-m+1.)*(l+m));
                V2[i] =-0.5*ipow(-(l-1))*sqrt((l+m+1.)*(l-m));
            }
        }


        COMPLEX rsnormal = stack->rs12(0.,lambda,vacuum,substrate);
        COMPLEX rpnormal = stack->rp12(0.,lambda,vacuum,substrate);

        vector<COMPLEX> reflect_s(npts);
        vector<COMPLEX> reflect_p(npts);
        vector<vector<COMPLEX> > Umatrix(npts,vector<COMPLEX>(index(LMAX,LMAX)+1));
        vector<vector<COMPLEX> > Vmatrix(npts,vector<COMPLEX>(index(LMAX,LMAX)+1));
        vector<vector<COMPLEX> > dminusmatrix(npts,vector<COMPLEX>(index(LMAX,LMAX)+1));
        vector<vector<COMPLEX> > dplusmatrix(npts,vector<COMPLEX>(index(LMAX,LMAX)+1));

        for (int j=0; j<npts; ++j) {
            ostringstream temp;
            temp << "Creating lookup tables: j = " << j << '/' << npts;
            message(temp.str());


            COMPLEX cosa = 1.0 - Gauss_Laguerre_Integration::zeros[ipt][j]/(2.*cI*qq);
            COMPLEX alpha = arccosine(cosa);
            COMPLEX cosa2 = cos((pi-alpha)/2.);
            COMPLEX sina2 = sin((pi-alpha)/2.);

            if (Norm_Inc_Approx) {
                reflect_s[j] = rsnormal;
                reflect_p[j] = rpnormal;
            } else {
				COMPLEX rs = stack->rs12(alpha, lambda, vacuum, substrate);
				COMPLEX rp = stack->rp12(alpha, lambda, vacuum, substrate);


                reflect_s[j] = (rs == rs) ? rs : 0.;
                reflect_p[j] = (rp == rp) ? rp : 0.;
            }


            for (int l=1; l<=LMAX; ++l) {
                for (int m=0; m<=l; ++m) {
                    int i = index(l,m);
                    Umatrix[j][i] = U1[i]*Ptilde(l-1,m-1,cosa)+
                                    U2[i]*Ptilde(l-1,m+1,cosa)+
                                    U3[i]*Ptilde(l+1,m-1,cosa)+
                                    U4[i]*Ptilde(l+1,m+1,cosa);
                    Vmatrix[j][i] = V1[i]*Ptilde(l,m-1,cosa)+
                                    V2[i]*Ptilde(l,m+1,cosa);

                    dminusmatrix[j][i] = dminus(l,m,cosa2,sina2);
                    dplusmatrix[j][i] = dplus(l,m,cosa2,sina2);
                }
            }
        }



        COMPLEX expx = exp(2.*cI*qq)/(2.*cI*qq);


        for (int m=0; m<=LMAX; ++m) {
            ostringstream temp;
            temp << "Constructing A matrix: m = " << m << '/' << LMAX;
            message(temp.str());

            int beginl = (m==0) ? 1 : m;
            for (int l=beginl; l<=LMAX; ++l) {
                for (int l_=beginl; l_<=LMAX; ++l_) {

                    COMPLEX Aee=0,Ahe=0,Aeh=0,Ahh=0;

                    int ii = index(l,m);
                    int ii_ = index(l_,m);


                    for (int _k=0; _k<npts; ++_k) {
                        double weight = Gauss_Laguerre_Integration::weights[ipt][_k];
                        Aee += weight*(reflect_p[_k]*Vmatrix[_k][ii]*dminusmatrix[_k][ii_]
                                       +reflect_s[_k]*Umatrix[_k][ii]*dplusmatrix[_k][ii_]);
                        Ahe += weight*(reflect_p[_k]*Vmatrix[_k][ii]*dplusmatrix[_k][ii_]
                                       +reflect_s[_k]*Umatrix[_k][ii]*dminusmatrix[_k][ii_]);
                        Aeh += weight*(reflect_p[_k]*Umatrix[_k][ii]*dminusmatrix[_k][ii_]
                                       +reflect_s[_k]*Vmatrix[_k][ii]*dplusmatrix[_k][ii_]);
                        Ahh += weight*(reflect_p[_k]*Umatrix[_k][ii]*dplusmatrix[_k][ii_]
                                       +reflect_s[_k]*Vmatrix[_k][ii]*dminusmatrix[_k][ii_]);
                    }


                    COMPLEX aa = expx*mpow(m-1)*ipow(l_-1)*lvector[l_];

                    Aee *= aa;
                    Ahe *= -cI*aa;
                    Aeh *= cI*aa;
                    Ahh *= aa;



                    _A[LMAX+m][2*(LMAX-l_)][2*(LMAX-l)] = Aee;
                    _A[LMAX+m][2*(LMAX-l_)+1][2*(LMAX-l)] = Ahe;
                    _A[LMAX+m][2*(LMAX-l_)+1][2*(LMAX-l)+1] = Ahh;
                    _A[LMAX+m][2*(LMAX-l_)][2*(LMAX-l)+1] = Aeh;




                    _A[LMAX-m][2*(LMAX-l_)][2*(LMAX-l)] = Aee;
                    _A[LMAX-m][2*(LMAX-l_)+1][2*(LMAX-l)] = -Ahe;
                    _A[LMAX-m][2*(LMAX-l_)+1][2*(LMAX-l)+1] = Ahh;
                    _A[LMAX-m][2*(LMAX-l_)][2*(LMAX-l)+1] = -Aeh;

                }
            }
        }
    }




#define VV(a,b) V[(a)][(b)]

    void
    BobVlieg_Supp::
    cvert(COMPLEX** V, int N)
    {
        vector<int> W(N);

        COMPLEX Y,Z;

        double S,T;
        int   I,J,K=0,L,M,P;

        if ( N != 1 ) {
            L = 0;
            M = 1;
Line10:
            if ( L != N ) {
                K = L;
                L = M;
                M = M + 1;



                P = L;
                S = abs(VV(L-1,L-1));
                if ( M <= N ) {
                    for (I=M; I<=N; ++I) {
                        T = abs(VV(I-1,L-1));
                        if ( T > S ) {
                            P = I;
                            S = T;
                        }
                    }
                    W[L-1] = P;
                }
                if ( S == 0. ) throw SCATMECH_exception("Matrix has no inverse (1)");
                Y = VV(P-1,L-1);
                VV(P-1,L-1) = VV(L-1,L-1);



                VV(L-1,L-1) = COMPLEX(-1.,0.);
                Y = COMPLEX(1.,0.)/Y;
                for (I=1; I<=N; ++I) {
                    VV(I-1,L-1) = -Y*VV(I-1,L-1);
                }
                J = L;

                while (1) {
                    do {
                        J = J + 1;
                        if ( J > N ) J = 1;
                        if ( J == L ) goto Line10;
                        Z = VV(P-1,J-1);
                        VV(P-1,J-1) = VV(L-1,J-1);
                        VV(L-1,J-1) = Z;
                    } while ( abs(Z) == 0. );



                    if ( K != 0 ) {
                        for (I=1; I<=K; ++I) {
                            VV(I-1,J-1) = VV(I-1,J-1) + Z*VV(I-1,L-1);
                        }
                    }
                    VV(L-1,J-1) = Y*Z;
                    if ( M <= N ) {
                        for (I=M; I<=N; ++I) {
                            VV(I-1,J-1) = VV(I-1,J-1) + Z*VV(I-1,L-1);
                        }
                    }
                };



            }

            do {
                L = W[K-1];
                for (I=1; I<=N; ++I) {
                    Y = VV(I-1,L-1);
                    VV(I-1,L-1) = VV(I-1,K-1);
                    VV(I-1,K-1) = Y;
                }
                K = K - 1;
            } while ( K > 0 );
            return;
        }
        if ( abs(VV(0,0)) != 0. ) {
            VV(0,0) = COMPLEX(1.,0.)/VV(0,0);
            return;
        }

        throw SCATMECH_exception("Matrix has no inverse (2)");
    }









    static
    void
    MatrixInvert_by_Series(COMPLEX **a,int n,int order)
    {
        vector<vector<COMPLEX> > temp(n,vector<COMPLEX>(n));
        vector<vector<COMPLEX> > prev(n,vector<COMPLEX>(n));
        vector<COMPLEX> diag(n);

        int i,j,k,q;


        for (i=0; i<n; ++i) {
            diag[i]=a[i][i];
        }


        for (i=0; i<n; ++i) {
            for (j=0; j<n; ++j) {
                a[i][j]/=diag[j];
            }
        }

        for (i=0; i<n; ++i) {
            for (j=0; j<n; ++j) {
                double unit_ij = ((i==j)?1.:0.);
                prev[i][j]=unit_ij;
                a[i][j]=unit_ij-a[i][j];
            }
        }

        for (q=0; q<order-1; ++q) {
            for (i=0; i<n; ++i) {
                for (j=0; j<n; ++j) {
                    double unit_ij=((i==j)?1.:0.);
                    temp[i][j]=0.;
                    for (k=0; k<n; ++k) {

                        temp[i][j]+= a[i][k]*prev[k][j];
                    }
                    temp[i][j]=unit_ij+temp[i][j];
                }
            }
            for (i=0; i<n; ++i) {
                for (j=0; j<n; ++j) {
                    prev[i][j]=temp[i][j];
                }
            }
        }

        for (i=0; i<n; ++i) {
            for (j=0; j<n; ++j) {
                a[i][j]=prev[i][j]/diag[i];
            }
        }
    }





    void
    Bobbert_Vlieger_BRDF_Model::
    invert_block_diagonal(ScatterTMatrix& AA)
    {
        int m;
        for (m=-LMAX; m<=LMAX; ++m) {
            ostringstream temp;
            temp << "Inverting matrix (m = " << m << ") ";
            message(temp.str());

            int n= 2*((m==0) ? LMAX : LMAX-abs(m)+1);

            if (order<0) {
                cvert(AA[LMAX+m],n);
            } else {
                MatrixInvert_by_Series(AA[LMAX+m],n,order);
            }
        }
        message("Done Inverting Matrices");
    }




    void
    Bobbert_Vlieger_BRDF_Model::
    setup()
    {
        Local_BRDF_Model::setup();

        N2 = sphere.index(lambda);

        int m,l,l_,f,f_;

        if (delta<0.)
            error("delta cannot be less than zero.");

        k   = 2.0*pi/lambda;
        q   = k*radius;
        qq  = k*(radius+delta);
        N0  = COMPLEX(substrate.index(lambda));

        r0 = radius;
        d = 0.;
        for (int i=0; i<spherecoat->get_n(); ++i) {
            r0 -= spherecoat->get_t()[i];
            d += spherecoat->get_t()[i];
        }
        if (r0<=0) throw SCATMECH_exception("Radius of particle, r, less than or equal to total thickness of coatings.");


        BH_LMAX = (int)(q+4.*pow(q,0.3333)+2.);


        if (lmax<=0) {
            LMAX = BH_LMAX-lmax;
        } else {
            LMAX = lmax;
        }

        if (BH_LMAX>LMAX) BH_LMAX=LMAX;

        sqrsize = index(LMAX,LMAX,1)+1;

        message("Allocating Memory");


        ScatMatrix.resize(LMAX);
        ScatMatrixInverse.resize(LMAX);
        A.resize(LMAX);

        Zp.resize(sqrsize);
        Zs.resize(sqrsize);
        eIP.resize(sqrsize);
        Wp.resize(sqrsize);
        Ws.resize(sqrsize);
        Vp.resize(sqrsize);
        Vs.resize(sqrsize);
        VSRp.resize(sqrsize);
        VSRs.resize(sqrsize);

        vector<COMPLEX> B(sqrsize);




        message("Calculating matrix A");
        Amatrix(A);


        message("Calculating matrix B");

        for (l=LMAX; l>=1; --l) {
            for (f=0; f<=1; ++f) {
                COMPLEX b = Bmatrix(l,0,f,l,0,f);
                for (m=-l; m<=l; ++m) {
                    B[index(l,m,f)]=b;
                }
            }
        }





        if (order!=0) {
            message("Calculating matrix (1/B)(1-BA)");
            for (m=-LMAX; m<=LMAX; ++m) {
                int beginl = (m==0) ? 1 : abs(m);
                for (l=beginl; l<=LMAX; ++l) {
                    for (f=0; f<=1; ++f) {
                        int i = index(l,m,f);
                        for (l_=beginl; l_<=LMAX; ++l_) {
                            for (f_=0; f_<=1; ++f_) {
                                int i_ = index(l_,m,f_);
                                int mm = LMAX+m;
                                int ll = 2*(LMAX-l)+f;
                                int ll_ = 2*(LMAX-l_)+f_;

                                ScatMatrix[mm][ll_][ll]=((double)(i==i_) -
                                                         B[i_]*A[mm][ll_][ll]);
                                ScatMatrix[mm][ll_][ll]/=B[i_];
                                ScatMatrixInverse[mm][ll_][ll]=ScatMatrix[mm][ll_][ll];
                            }
                        }
                    }
                }
            }




            invert_block_diagonal(ScatMatrix);
        } else {
            message("Calculating matrix (1/B)(1-BA)");
            for (m=-LMAX; m<=LMAX; ++m) {
                int beginl = (m==0) ? 1 : abs(m);
                for (l=beginl; l<=LMAX; ++l) {
                    for (f=0; f<=1; ++f) {
                        int i = index(l,m,f);
                        for (l_=beginl; l_<=LMAX; ++l_) {
                            for (f_=0; f_<=1; ++f_) {
                                int i_ = index(l_,m,f_);
                                int mm = LMAX+m;
                                int ll = 2*(LMAX-l)+f;
                                int ll_ = 2*(LMAX-l_)+f_;

                                ScatMatrix[mm][ll_][ll] = (i==i_ ? B[i_] : 0.);
                                ScatMatrixInverse[mm][ll_][ll] = (i==i_ ? 1./B[i_] : 0.);
                            }
                        }
                    }
                }
            }
        }

        old_thetai=-1000;
        old_thetas=-1000;
        old_phis=-1000;

    }





    void
    Bobbert_Vlieger_BRDF_Model::
    set_geometry(double _thetai,double _thetas,double _phis)
    {
        SETUP();

        if (abs(_thetai)<1E-6) _thetai=1E-6;
        if (abs(_thetas)<1E-6) _thetas=1E-6;

        if (_thetai<0) {
            _thetai=-_thetai;
            _phis = pi+_phis;
        }

        if (_thetas<0) {
            _thetas=-_thetas;
            _phis=pi+_phis;
        }

        if (_thetai!=old_thetai) {
            calculate_W(_thetai);
            old_thetai=_thetai;
        }

        if (_thetas!=old_thetas) {
            calculate_Z(pi-_thetas);
            old_thetas=_thetas;
        }

        if (_phis!=old_phis) {
            calculate_eIP(_phis);
            old_phis=_phis;
        }
    }



    COMPLEX
    Bobbert_Vlieger_BRDF_Model::
    Epp(double _thetai,double _thetas,double _phis)
    {
        set_geometry(_thetai,_thetas,_phis);
        return E(Wp,Zp);
    }



    COMPLEX
    Bobbert_Vlieger_BRDF_Model::
    Eps(double _thetai,double _thetas,double _phis)
    {
        set_geometry(_thetai,_thetas,_phis);
        return E(Wp,Zs);
    }



    COMPLEX
    Bobbert_Vlieger_BRDF_Model::
    Esp(double _thetai,double _thetas,double _phis)
    {
        set_geometry(_thetai,_thetas,_phis);
        return E(Ws,Zp);
    }



    COMPLEX
    Bobbert_Vlieger_BRDF_Model::
    Ess(double _thetai,double _thetas,double _phis)
    {
        set_geometry(_thetai,_thetas,_phis);
        return E(Ws,Zs);
    }




    CVector
    Bobbert_Vlieger_BRDF_Model::
    EField(double _thetai,const JonesVector& inpol,const Vector R)
    {
        double rr = sqrt(sqr(R.x)+sqr(R.y)+sqr(R.z));
        double theta = acos(R.z/rr);
        double phi = atan2(R.y,R.x);

        if (rr<radius || R.z > radius+delta) return CVector(0,0,0);
        if (rr>radius+delta) return CVector(0,0,0);

        set_geometry(_thetai,0.,0.);

        COMPLEX Esr = 0;
        COMPLEX Estheta = 0;
        COMPLEX Esphi = 0;
        COMPLEX Epr = 0;
        COMPLEX Eptheta = 0;
        COMPLEX Epphi = 0;

        double costheta = cos(theta);
        double sintheta = sin(theta);
        double tantheta = tan(theta);

        for (int i=0,l=1; l<=LMAX; ++l) {
            for (int m=-l; m<=l; ++m) {
                COMPLEX wse = Ws[i];
                COMPLEX wpe = Wp[i];
                COMPLEX vse = VSRs[i]+Vs[i];
                COMPLEX vpe = VSRp[i]+Vp[i];
                i++;
                COMPLEX wsh = Ws[i];
                COMPLEX wph = Wp[i];
                COMPLEX vsh = VSRs[i]+Vs[i];
                COMPLEX vph = VSRp[i]+Vp[i];
                i++;

                COMPLEX eImphi = exp(COMPLEX(0,m)*phi);
                COMPLEX Llm = Legendre(l,m,costheta);
                COMPLEX Llmp1 = Legendre(l,m+1,costheta);
                COMPLEX hlkr = h(l,k*rr);
                COMPLEX hlp1kr = h(l+1,k*rr);
                COMPLEX jlkr = j(l,k*rr);
                COMPLEX jlp1kr =j(l+1,k*rr);
                double sqrt12l = sqrt(1.+2.*l);

                COMPLEX ere = ((double)l*(1. + l)*eImphi*hlkr*Llm*mpow(2*m)*sqrt12l*sqrtFact(l - m))/
                              (rr*sqrtFact(l + m));
                COMPLEX erh = 0;
                COMPLEX ethetae = (eImphi*((1. + l)*hlkr - k*rr*hlp1kr)*mpow(2*m)*sqrt12l*
                                   sqrtFact(l - m)*(-Llmp1 + ((double)m*Llm)/tantheta))/
                                  (rr*sqrtFact(l + m));
                COMPLEX ethetah = -((k*(double)m*eImphi*hlkr*Llm*mpow(2*m)*sqrt12l*
                                     sqrtFact(l - m))/(sqrtFact(l + m)*sintheta));
                COMPLEX ephie = (COMPLEX(0,1)*(double)m*eImphi*((1. + l)*hlkr - k*rr*hlp1kr)*
                                 Llm*mpow(2*m)*sqrt12l*sqrtFact(l - m))/(rr*sqrtFact(l + m)*sintheta);
                COMPLEX ephih = (COMPLEX(0,1)*k*eImphi*hlkr*mpow(2*m)*sqrt12l*sqrtFact(l - m)*
                                 (Llmp1 - ((double)m*Llm)/tantheta))/sqrtFact(l + m);

                Esr += ( ere*wse + erh*wsh );
                Epr += ( ere*wpe + erh*wph );
                Estheta += ( ethetae*wse + ethetah*wsh );
                Eptheta += ( ethetae*wpe + ethetah*wph );
                Esphi += ( ephie*wse + ephih*wsh );
                Epphi += ( ephie*wpe + ephih*wph );

                ere = ((double)l*(1. + l)*eImphi*jlkr*Llm*mpow(2*m)*sqrt12l*sqrtFact(l - m))/
                      (rr*sqrtFact(l + m));
                erh = 0;
                ethetae = (eImphi*((1. + l)*jlkr - k*rr*jlp1kr)*mpow(2*m)*sqrt12l*
                           sqrtFact(l - m)*(-Llmp1 + ((double)m*Llm)/tantheta))/
                          (rr*sqrtFact(l + m));
                ethetah = -((k*(double)m*eImphi*jlkr*Llm*mpow(2*m)*sqrt12l*
                             sqrtFact(l - m))/(sqrtFact(l + m)*sintheta));
                ephie = (COMPLEX(0,1)*(double)m*eImphi*((1. + l)*jlkr - k*rr*jlp1kr)*
                         Llm*mpow(2*m)*sqrt12l*sqrtFact(l - m))/(rr*sqrtFact(l + m)*sintheta);
                ephih = (COMPLEX(0,1)*k*eImphi*jlkr*mpow(2*m)*sqrt12l*sqrtFact(l - m)*
                         (Llmp1 - ((double)m*Llm)/tantheta))/sqrtFact(l + m);

                Esr += ( ere*vse + erh*vsh );
                Epr += ( ere*vpe + erh*vph );
                Estheta += ( ethetae*vse + ethetah*vsh );
                Eptheta += ( ethetae*vpe + ethetah*vph );
                Esphi += ( ephie*vse + ephih*vsh );
                Epphi += ( ephie*vpe + ephih*vph );
            }
        }

        COMPLEX Er = Esr*inpol.S() + Epr*inpol.P();
        COMPLEX Etheta = Estheta*inpol.S() + Eptheta*inpol.P();
        COMPLEX Ephi = Esphi*inpol.S() + Epphi*inpol.P();

        return Er*CVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)) +
               Etheta*CVector(cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)) +
               Ephi*CVector(-sin(phi),cos(phi),0);
    }




    JonesMatrix
    Bobbert_Vlieger_BRDF_Model::
    jonesDSC()
    {
        set_geometry(thetai,thetas,phis);

        COMPLEX pp=E(Wp,Zp);
        COMPLEX ps=E(Wp,Zs);
        COMPLEX sp=E(Ws,Zp);
        COMPLEX ss=E(Ws,Zs);

        return JonesMatrix(pp,ss,ps,sp);
    }




    Bobbert_Vlieger_BRDF_Model::
    Bobbert_Vlieger_BRDF_Model()
    {

        LMAX=0;


        ScatMatrix=0;
        ScatMatrixInverse=0;


        old_thetai=-1000.;
        old_thetas=-1000.;
        old_phis=-1000;
    }

    MuellerMatrix
    Bobbert_Vlieger_BRDF_Model::
    Specular(double theta)
    {
        SETUP();







        switch (type) {
        case 0:
        {
            JonesMatrix X = JonesDSC(theta, theta, 0., 0.);
            JonesMatrix r = stack->r12(theta, lambda, vacuum, substrate);
            r *= exp(2. * cI * qq * cos(theta));
            MuellerMatrix sigma = (4. * pi / k) * ReCrossMueller(X, r);

            return MuellerMatrix(r) - sigma * (density / cos(theta));
        }
        break;
        case 1:
        {
            double index = substrate.n(lambda);
            double sint = sin(theta) / index;
            if (sint >= 1. || substrate.k(lambda) != 0) return MuellerZero();
            double thetat = asin(sint);
            JonesMatrix X = JonesDSC(theta, thetat, 0., 0.);
            COMPLEX phase = exp(cI * qq * (cos(theta) - index * cos(thetat)));

            JonesMatrix t = stack->t12(theta, lambda, vacuum, substrate);
            t *= phase;
            MuellerMatrix sigma = (4. * pi / k / sqrt(index)) * ReCrossMueller(X, t);

            return MuellerMatrix(t) * (cos(thetat) / cos(theta) * index) - sigma * (density / cos(theta));
        }
        break;
        case 2:
        {
            double index = substrate.n(lambda);
            JonesMatrix X = JonesDSC(theta, theta, 0., 0.);
            JonesMatrix r = stack->r21i(theta, lambda, substrate, vacuum);
            r *= exp(-2. * cI * index * qq * cos(theta));
            MuellerMatrix sigma = (4. * pi / k / index) * ReCrossMueller(X, r);

            return MuellerMatrix(r) - sigma * (density / cos(theta));
        }
        break;
        case 3:
        {
            double index = substrate.n(lambda);
            double sint = sin(theta) * index;
            COMPLEX cost = sqrt(1. - sqr(sint));
            if (imag(cost) < 0) cost = -cost;
            if (sint >= 1.) return MuellerZero();
            double thetat = asin(sint);
            JonesMatrix X = JonesDSC(theta, thetat, 0., 0.);
            COMPLEX phase = exp(cI * qq * (cost - index * cos(theta)));

            JonesMatrix t = stack->t21i(theta, lambda, substrate, vacuum);
            t *= phase;
            MuellerMatrix sigma = (4. * pi / k / sqrt(index)) * ReCrossMueller(X, t);

            return MuellerMatrix(t) * (cos(theta) / index / cos(thetat)) - sigma * (density / cos(theta));
        }
        break;
        default:
            error("Invalid type = " + to_string(type));
        }
        return MuellerZero();
    }

    DEFINE_MODEL(Bobbert_Vlieger_BRDF_Model,Local_BRDF_Model,
                 "Theory for scattering by a sphere on a substrate.");

    DEFINE_PARAMETER(Bobbert_Vlieger_BRDF_Model,dielectric_function,sphere,"Sphere optical properties","(1.59,0)",0xFF);
    DEFINE_PARAMETER(Bobbert_Vlieger_BRDF_Model,double,radius,"Particle radius [um]","0.05",0xFF);
    DEFINE_PTRPARAMETER(Bobbert_Vlieger_BRDF_Model,StackModel_Ptr,spherecoat,"Coatings on the sphere","No_StackModel",0xFF);
    DEFINE_PTRPARAMETER(Bobbert_Vlieger_BRDF_Model,StackModel_Ptr,stack,"Substrate films","No_StackModel",0xFF);
    DEFINE_PARAMETER(Bobbert_Vlieger_BRDF_Model,double,delta,"Separation of particle from substrate [um] (in contact: 0)","0",0xFF);
    DEFINE_PARAMETER(Bobbert_Vlieger_BRDF_Model,int,lmax,"Maximum spherical harmonic order (use Bohren & Huffman estimate: 0)","0",0xFF);
    DEFINE_PARAMETER(Bobbert_Vlieger_BRDF_Model,int,order,"Perturbation order (exact: -1)","-1",0xFF);
    DEFINE_PARAMETER(Bobbert_Vlieger_BRDF_Model,int,Norm_Inc_Approx,"Normal Incidence Approximation (exact: 0)","0",0xFF);
    DEFINE_PARAMETER(Bobbert_Vlieger_BRDF_Model,int,improve,"Iterative improvement steps (recommend: 3)","3",0xFF);

}



