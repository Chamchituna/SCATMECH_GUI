











#include <vector>
#include <cmath>
#include <limits>
#include "scatmech.h"

#include "matrixmath.h"

using namespace std;

#define THROW(x) throw SCATMECH_exception(x);

namespace SCATMECH {

    namespace CMLIB {









        inline double CABS1(COMPLEX ZDUM) {
            return abs(real(ZDUM)) + abs(imag(ZDUM));
        }
        inline COMPLEX CSIGN(const COMPLEX& ZDUM1,const COMPLEX& ZDUM2) {
            return abs(ZDUM1)*(ZDUM2/abs(ZDUM2));
        }
        inline int MIN0(int a,int b) {
            return a>b? b : a;
        }
        inline int MAX0(int a,int b) {
            return a>b? a : b;
        }
        template<class T> T max(const T& a,const T& b) {
            return a>b ? a : b;
        }
        inline double AMAX1(double a,double b,double c,double d,double e) {
            return max(a,max(b,max(c,max(d,e))));
        }


        void CSVDC(CFARRAY X,const int& LDX,const int& N,const int& P,CFARRAY S,CFARRAY E,CFARRAY U,const int& LDU,CFARRAY V,const int& LDV,CFARRAY WORK,const int& JOB,int& INFO)
        {

            X.array(LDX,P);
            if (JOB/10==1) U.array(LDU,N);
            else if (JOB/10==2) U.array(LDU,N<P?N:P);
            else U.array(LDU);
            V.array(LDV,P);





















































































































            using namespace SCATMECH::CMLIB;
            int I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,
                MM,MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1;
            COMPLEX T,R;
            double B,C,CS,EL,EMM1,F,G,SCALE,SHIFT,SL,SM,SN,SMM1,T1,TEST,ZTEST;
            bool WANTU,WANTV;
            COMPLEX ZDUM,ZDUM1,ZDUM2;



	    int MAXIT = 30;





            WANTU = false;
            WANTV = false;
            JOBU = (JOB%100)/10;
            NCU = N;
            if (JOBU > 1) NCU = MIN0(N,P);
            if (JOBU != 0) WANTU = true;
            if ((JOB%10) != 0) WANTV = true;




            INFO = 0;
            NCT = MIN0(N-1,P);
            NRT = MAX0(0,MIN0(P-2,N));
            LU = MAX0(NCT,NRT);
            if (LU >= 1) {
                for (L=1; L<=LU; ++L) {
                    LP1 = L + 1;
                    if (L <= NCT) {




                        S(L) = COMPLEX(SCNRM2(N-L+1,X(L,L),1),0.0E0);
                        if (CABS1(S(L)) != 0.0E0) {
                            if (CABS1(X(L,L)) != 0.0E0) S(L) = CSIGN(S(L),X(L,L));
                            CSCAL(N-L+1,1.0E0/S(L),X(L,L),1);
                            X(L,L) = COMPLEX(1.0E0,0.0E0) + X(L,L);
                        }
                        S(L) = -S(L);
                    }
                    if (P >= LP1) {
                        for (J = LP1; J<=P; ++J) {
                            if (L <= NCT) {
                                if (CABS1(S(L)) != 0.0E0) {



                                    T = -CDOTC(N-L+1,X(L,L),1,X(L,J),1)/X(L,L);
                                    CAXPY(N-L+1,T,X(L,L),1,X(L,J),1);
                                }
                            }



                            E(J) = conj(X(L,J));
                        }
                    }
                    if (!(!WANTU || L > NCT)) {




                        for (I=L; I<=N; ++I) {
                            U(I,L) = X(I,L);
                        }
                    }
                    if (L <= NRT) {




                        E(L) = COMPLEX(SCNRM2(P-L,E(LP1),1),0.0E0);
                        if (CABS1(E(L)) != 0.0E0) {
                            if (CABS1(E(LP1)) != 0.0E0) E(L) = CSIGN(E(L),E(LP1));
                            CSCAL(P-L,1.0E0/E(L),E(LP1),1);
                            E(LP1) = COMPLEX(1.0E0,0.0E0) + E(LP1);
                        }
                        E(L) = -conj(E(L));
                        if (!(LP1 > N || CABS1(E(L)) == 0.0E0)) {



                            for (I=LP1; I<=N; ++I) {
                                WORK(I) = COMPLEX(0.0E0,0.0E0);
                            }
                            for (J=LP1; J<=P; ++J) {
                                CAXPY(N-L,E(J),X(LP1,J),1,WORK(LP1),1);
                            }
                            for (J=LP1; J<=P; ++J) {
                                COMPLEX temp = conj(-E(J)/E(LP1));
                                CAXPY(N-L,temp,WORK(LP1),1,X(LP1,J),1);
                            }
                        }
                        if (WANTV) {




                            for (I=LP1; I<=P; ++I) {
                                V(I,L) = E(I);
                            }
                        }
                    }
                }
            }



            M = MIN0(P,N+1);
            NCTP1 = NCT + 1;
            NRTP1 = NRT + 1;
            if (NCT < P) S(NCTP1) = X(NCTP1,NCTP1);
            if (N < M) S(M) = COMPLEX(0.0E0,0.0E0);
            if (NRTP1 < M) E(NRTP1) = X(NRTP1,M);
            E(M) = COMPLEX(0.0E0,0.0E0);



            if (WANTU) {
                if (NCU >= NCTP1) {
                    for (J=NCTP1; J<=NCU; ++J) {
                        for (I=1; I<=N; ++I) {
                            U(I,J) = COMPLEX(0.0E0,0.0E0);
                        }
                        U(J,J) = COMPLEX(1.0E0,0.0E0);
                    }
                }
                if (NCT >= 1) {
                    for (LL = 1; LL<=NCT; ++LL) {
                        L = NCT - LL + 1;
                        if (CABS1(S(L)) != 0.0E0) {
                            LP1 = L + 1;
                            if (NCU >= LP1) {
                                for (J=LP1; J<=NCU; ++J) {
                                    T = -CDOTC(N-L+1,U(L,L),1,U(L,J),1)/U(L,L);
                                    CAXPY(N-L+1,T,U(L,L),1,U(L,J),1);
                                }
                            }
                            CSCAL(N-L+1,COMPLEX(-1.0E0,0.0E0),U(L,L),1);
                            U(L,L) = COMPLEX(1.0E0,0.0E0) + U(L,L);
                            LM1 = L - 1;
                            if (LM1 >= 1) {
                                for (I = 1; I<=LM1; ++I) {
                                    U(I,L) = COMPLEX(0.0E0,0.0E0);
                                }
                            }
                        } else {
                            for (I=1; I<=N; ++I) {
                                U(I,L) = COMPLEX(0.0E0,0.0E0);
                            }
                            U(L,L) = COMPLEX(1.0E0,0.0E0);
                        }
                    }
                }
            }



            if (WANTV) {
                for (LL=1; LL<=P; ++LL) {
                    L = P - LL + 1;
                    LP1 = L + 1;
                    if (L <= NRT) {
                        if (CABS1(E(L)) != 0.0E0) {
                            for (J = LP1; J<=P; ++J) {
                                T = -CDOTC(P-L,V(LP1,L),1,V(LP1,J),1)/V(LP1,L);
                                CAXPY(P-L,T,V(LP1,L),1,V(LP1,J),1);
                            }
                        }
                    }
                    for (I=1; I<=P; ++I) {
                        V(I,L) = COMPLEX(0.0E0,0.0E0);
                    }
                    V(L,L) = COMPLEX(1.0E0,0.0E0);
                }
            }



            for (I=1; I<=M; ++I) {
                if (CABS1(S(I)) != 0.0E0) {
                    T = COMPLEX(abs(S(I)),0.0E0);
                    R = S(I)/T;
                    S(I) = T;
                    if (I < M) E(I) = E(I)/R;
                    if (WANTU) CSCAL(N,R,U(1,I),1);
                }

                if (I == M) goto Line390;
                if (CABS1(E(I)) != 0.0E0) {
                    T = COMPLEX(abs(E(I)),0.0E0);
                    R = T/E(I);
                    E(I) = T;
                    S(I+1) = S(I+1)*R;
                    if (WANTV) CSCAL(P,R,V(1,I+1),1);
                }
            }

Line390:
            ;



            MM = M;
            ITER = 0;



            while (M != 0) {




                if (ITER >= MAXIT) {
                    INFO = M;

                    return;
                }











                for (LL=1; LL<=M; ++LL) {
                    L = M - LL;

                    if (L == 0) goto Line440;
                    TEST = abs(S(L)) + abs(S(L+1));
                    ZTEST = TEST + abs(E(L));
                    if (ZTEST == TEST) {
                        E(L) = COMPLEX(0.0E0,0.0E0);

                        goto Line440;
                    }
                }
Line440:
                if (L == M - 1) {
                    KASE = 4;
                } else {
                    LP1 = L + 1;
                    MP1 = M + 1;
                    for (LLS = LP1; LLS<=MP1; ++LLS) {
                        LS = M - LLS + LP1;

                        if (LS == L) goto Line480;
                        TEST = 0.0E0;
                        if (LS != M) TEST = TEST + abs(E(LS));
                        if (LS != L + 1) TEST = TEST + abs(E(LS-1));
                        ZTEST = TEST + abs(S(LS));
                        if (ZTEST == TEST) {
                            S(LS) = COMPLEX(0.0E0,0.0E0);

                            goto Line480;
                        }
                    }
Line480:
                    if (LS == L) {
                        KASE = 3;
                    } else if (LS == M) {
                        KASE = 1;
                    } else {
                        KASE = 2;
                        L = LS;
                    }
                }
                L = L + 1;



                switch (KASE) {



                    case 1:
                        MM1 = M - 1;
                        F = real(E(M-1));
                        E(M-1) = COMPLEX(0.0E0,0.0E0);
                        for (KK=L; KK<=MM1; ++KK) {
                            K = MM1 - KK + L;
                            T1 = real(S(K));
                            SROTG(T1,F,CS,SN);
                            S(K) = COMPLEX(T1,0.0E0);
                            if (K != L) {
                                F = -SN*real(E(K-1));
                                E(K-1) = CS*E(K-1);
                            }
                            if (WANTV) CSROT(P,V(1,K),1,V(1,M),1,CS,SN);
                        }
                        break;



                    case 2:
                        F = real(E(L-1));
                        E(L-1) = COMPLEX(0.0E0,0.0E0);
                        for (K=L; K<=M; ++K) {
                            T1 = real(S(K));
                            SROTG(T1,F,CS,SN);
                            S(K) = COMPLEX(T1,0.0E0);
                            F = -SN*real(E(K));
                            E(K) = CS*E(K);
                            if (WANTU) CSROT(N,U(1,K),1,U(1,L-1),1,CS,SN);
                        }
                        break;



                    case 3:



                        SCALE = AMAX1(abs(S(M)),abs(S(M-1)),abs(E(M-1)),abs(S(L)),abs(E(L)));
                        SM = real(S(M))/SCALE;
                        SMM1 = real(S(M-1))/SCALE;
                        EMM1 = real(E(M-1))/SCALE;
                        SL = real(S(L))/SCALE;
                        EL = real(E(L))/SCALE;
                        B = ((SMM1 + SM)*(SMM1 - SM) + sqr(EMM1))/2.0E0;
                        C = sqr(SM*EMM1);
                        SHIFT = 0.0E0;
                        if (!(B == 0.0E0 && C == 0.0E0)) {
                            SHIFT = sqrt(sqr(B)+C);
                            if (B < 0.0E0) SHIFT = -SHIFT;
                            SHIFT = C/(B + SHIFT);
                        }
                        F = (SL + SM)*(SL - SM) - SHIFT;
                        G = SL*EL;



                        MM1 = M - 1;
                        for (K=L; K<=MM1; ++K) {
                            SROTG(F,G,CS,SN);
                            if (K != L) E(K-1) = COMPLEX(F,0.0E0);
                            F = CS*real(S(K)) + SN*real(E(K));
                            E(K) = CS*E(K) - SN*S(K);
                            G = SN*real(S(K+1));
                            S(K+1) = CS*S(K+1);
                            if (WANTV) CSROT(P,V(1,K),1,V(1,K+1),1,CS,SN);
                            SROTG(F,G,CS,SN);
                            S(K) = COMPLEX(F,0.0E0);
                            F = CS*real(E(K)) + SN*real(S(K+1));
                            S(K+1) = -SN*E(K) + CS*S(K+1);
                            G = SN*real(E(K+1));
                            E(K+1) = CS*E(K+1);
                            if (WANTU && K < N) CSROT(N,U(1,K),1,U(1,K+1),1,CS,SN);
                        }
                        E(M-1) = COMPLEX(F,0.0E0);
                        ITER = ITER + 1;
                        break;



                    case 4:



                        if (real(S(L)) < 0.0E0) {
                            S(L) = -S(L);
                            if (WANTV) CSCAL(P,COMPLEX(-1.0E0,0.0E0),V(1,L),1);
                        }




                        while (L != MM) {

                            if (real(S(L)) >= real(S(L+1))) goto Line640;
                            T = S(L);
                            S(L) = S(L+1);
                            S(L+1) = T;
                            if (WANTV && L < P) CSWAP(P,V(1,L),1,V(1,L+1),1);
                            if (WANTU && L < N) CSWAP(N,U(1,L),1,U(1,L+1),1);
                            L = L + 1;
                        }
Line640:
                        ITER = 0;
                        M = M - 1;
                        break;
                }
            }
            return;
        }

        void CSROT(const int& N,CFARRAY CX,const int& INCX,CFARRAY CY,const int& INCY,const double& C,const double& S)
        {
            COMPLEX CTEMP;
            int I,IX,IY;












































            if (N<=0) return;
            if (INCX==1&&INCY==1) goto Line20;




            IX = 1;
            IY = 1;
            if (INCX<0) IX = (-N+1)*INCX + 1;
            if (INCY<0) IY = (-N+1)*INCY + 1;

            for (I=1; I<=N; ++I) {
                CTEMP = C*CX(IX) + S*CY(IY);
                CY(IY) = C*CY(IY) - S*CX(IX);
                CX(IX) = CTEMP;
                IX = IX + INCX;
                IY = IY + INCY;
            }
            return;



Line20:
            for (I=1; I<=N; ++I) {
                CTEMP = C*CX(I) + S*CY(I);
                CY(I) = C*CY(I) - S*CX(I);
                CX(I) = CTEMP;
            }
            return;
        }

        double SCNRM2(const int& N,CFARRAY CX,const int& INCX)
        {





































































            bool IMAG, SCALE;
            int          NEXT,NN,I;
            double       HITEST, SUM, XMAX, ABSX, _SCNRM2;
            const double ONE = 1.;
            const double ZERO = 0.;
            const double CUTLO = sqrt(std::numeric_limits<double>::min()/std::numeric_limits<double>::epsilon());
            const double CUTHI = sqrt(std::numeric_limits<double>::max());

            if (N > 0) goto Line10;
            _SCNRM2 = ZERO;
            goto Line300;

Line10:
            NEXT = 30;
            SUM = ZERO;
            NN = N * INCX;

            for (I=1; I<=NN; I+=INCX) {
                ABSX = fabs(real(CX(I)));
                IMAG = false;
                if (NEXT==30) goto Line30;
                if (NEXT==50) goto Line50;
                if (NEXT==70) goto Line70;
                if (NEXT==90) goto Line90;
                if (NEXT==110) goto Line110;
                throw;
Line30:
                if ( ABSX > CUTLO) goto Line85;
                NEXT = 50;
                SCALE = false;



Line50:
                if ( ABSX == ZERO) goto Line200;
                if ( ABSX > CUTLO) goto Line85;


                NEXT = 70;
                goto Line105;



Line100:
                NEXT = 110;
                SUM = (SUM / ABSX) / ABSX;
Line105:
                SCALE = true;
                XMAX = ABSX;
                goto Line115;




Line70:
                if ( ABSX > CUTLO ) goto Line75;




Line110:
                if ( ABSX <= XMAX ) goto Line115;
                SUM = ONE + SUM * sqr(XMAX / ABSX);
                XMAX = ABSX;
                goto Line200;

Line115:
                SUM = SUM + sqr(ABSX/XMAX);
                goto Line200;




Line75:
                SUM = (SUM * XMAX) * XMAX;

Line85:
                NEXT = 90;
                SCALE = false;




                HITEST = CUTHI/double( N );



Line90:
                if (ABSX >= HITEST) goto Line100;
                SUM = SUM + sqr(ABSX);
Line200:
                ;


                if (IMAG) goto Line210;
                ABSX = fabs(imag(CX(I)));
                IMAG = true;
                if (NEXT==50) goto Line50;
                if (NEXT==70) goto Line70;
                if (NEXT==90) goto Line90;
                if (NEXT==110) goto Line110;
                throw;
Line210:
                ;
            }




            _SCNRM2 = sqrt(SUM);
            if (SCALE) _SCNRM2 = _SCNRM2 * XMAX;
Line300:
            return (_SCNRM2);
        }

        void SROTG(double& SA,double& SB,double& SC,double& SS)
        {






















































            double U,V,R;
            if (fabs(SA) <= fabs(SB)) goto Line10;



            U = SA + SA;
            V = SB / U;



            R = sqrt(.25 + sqr(V)) * U;



            SC = SA / R;
            SS = V * (SC + SC);
            SB = SS;
            SA = R;
            return;



Line10:
            if (SB == 0.) goto Line20;
            U = SB + SB;
            V = SA / U;




            SA = sqrt(.25 + sqr(V)) * U;



            SS = SB / SA;
            SC = V * (SS + SS);
            if (SC == 0.) goto Line15;
            SB = 1. / SC;
            return;
Line15:
            SB = 1.;
            return;



Line20:
            SC = 1.;
            SS = 0.;
            return;

        }
    }
}

