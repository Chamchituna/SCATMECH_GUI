











#include <vector>
#include <cmath>
#include "scatmech.h"

#include "matrixmath.h"




#ifdef USING_BLAS
extern "C" {
    struct cmplx {
        double re;
        double im;
    };
    extern void dcopy_(int *N,double *SX,int *INCX,double *SY,int *INCY);
    extern void zaxpy_(int *N, cmplx *CA, cmplx *CX, int *INCX, cmplx *CY, int *INCY);
    extern void zscal_(int *N, cmplx *CA, cmplx *CX, int *INCX);
    extern int izamax_(int *N, cmplx *CX, int *INCX);
    extern cmplx zdotc_(int *N, cmplx *CX, int *INCX, cmplx *CY, int *INCY);
    extern void zswap_(int *N, cmplx *CX, int *INCX, cmplx *CY, int *INCY);
    extern void zgemv_(const char *TRANS,int* M, int* N, cmplx *ALPHA, cmplx *A, int* LDA, cmplx *X, int* INCX, cmplx* BETA, cmplx* Y, int* INCY);
    extern void zgemm_(const char *TRANSA,const char*  TRANSB, int* M, int* N, int* K, cmplx* ALPHA, cmplx* A, int* LDA, cmplx* B, int* LDB, cmplx* BEAT, cmplx* C, int *LDC);

    extern void dgefa_(double* A,int* LDA,int* N,int* IPVT,int* INFO);
    extern void daxpy_(int* N,double* DA,double* DX,int* INCX,double* DY,int* INCY);
    extern void dscal_(int* N,double* DA,double* DX,int* INCX);
    extern int idamax_(int* N,double* DX,int* INCX);
    extern void dswap_(int* N,double* DX,int* INCX,double* DY,int* INCY);
    extern void dgedi_(double* A,int* LDA,int* N,int* IPVT,double* DET,double* WORK,int* JOB);
    extern void dgesl_(double& A,int* LDA,int* N,int* IPVT,double* B,int* JOB);
    extern double ddot_(int* N,double* DX,int* INCX,double* DY,int* INCY);
}
#endif

using namespace std;

#define THROW(x) throw SCATMECH_exception(x);

namespace SCATMECH {


    int eigen(vector<COMPLEX>& A,
              vector<COMPLEX>& Q,
              vector<COMPLEX>& W,int nmat)
    {






        return eigen(CFARRAY(&A[0]),CFARRAY(&Q[0]),CFARRAY(&W[0]),nmat);
    }

    void LUdecompose(std::vector<COMPLEX>& matrix,int nmat,std::vector<int>& pivot)
    {
        int info;
        CMLIB::CGEFA(matrix,nmat,nmat,pivot,info);
        if (info!=0) THROW("Singular matrix in LUdecompose");
    }

    void LUbacksubstitute(std::vector<COMPLEX>& matrix,int nmat,std::vector<int>& pivot,std::vector<COMPLEX>& b)
    {
        CMLIB::CGESL(matrix,nmat,nmat,pivot,b,0);
    }

    void LUImprove(vector<COMPLEX>& a, vector<COMPLEX>& alud, int n, vector<int>& indx, vector<COMPLEX>& b, vector<COMPLEX>& x)
    {
        vector<COMPLEX> r(n);

        for (int i=0; i<n; ++i) {
            r[i] = -b[i];
            for (int j=0; j<n; ++j) r[i] += a[i+n*j]*x[j];
        }
        LUbacksubstitute(alud,n,indx,r);
        for (int ii=1; ii<n; ++ii) x[ii] -= r[ii];
    }


    void Inverse(std::vector<COMPLEX>& matrix, int nmat)
    {
        using namespace CMLIB;

        CFARRAY work;
        work.allocate(nmat);
        IFARRAY pivot;
        pivot.allocate(nmat);
        CFARRAY det;
        det.allocate(2);
        int info;
        CGEFA(matrix,nmat,nmat,pivot,info);
        if (info!=0) THROW("Singular matrix in LUdecompose");
        CGEDI(matrix,nmat,nmat,pivot,det,work,1);
    }

    int eigen(CFARRAY A, CFARRAY Q, CFARRAY W, int nmat) {
        int info;
        DFARRAY work(3*nmat,1);
        CMLIB::CGEEV((double*)&(A[0]),nmat,nmat,(double*)&(Q[0]),(double*)&(W[0]),nmat,work,1,info);
        if (info!=0) THROW("Eigenvalue solver failed to converge");
        return info;
    }

    void LUdecompose(CFARRAY matrix,int nmat,IFARRAY pivot)
    {
        int info;
        CMLIB::CGEFA(matrix,nmat,nmat,pivot,info);
        if (info!=0) THROW("Singular matrix in LUdecompose");
    }

    void LUbacksubstitute(CFARRAY matrix,int nmat,IFARRAY pivot,CFARRAY b)
    {
        CMLIB::CGESL(matrix,nmat,nmat,pivot,b,0);
    }

    void LUdecompose(DFARRAY matrix,int nmat,IFARRAY pivot)
    {
        int info;
        CMLIB::DGEFA(matrix,nmat,nmat,pivot,info);
        if (info!=0) THROW("Singular matrix in LUdecompose");
    }

    void LUbacksubstitute(DFARRAY matrix,int nmat,IFARRAY pivot,DFARRAY b)
    {
        CMLIB::DGESL(matrix,nmat,nmat,pivot,b,0);
    }

    void LUImprove(CFARRAY a, CFARRAY alud, int n, IFARRAY indx, CFARRAY b, CFARRAY x)
    {
        CFARRAY r(n,1);

        for (int i=1; i<=n; ++i) {
            r(i) = -b(i);
            for (int j=1; j<=n; ++j) r(i) += a(i+n*j)*x(j);
        }
        LUbacksubstitute(alud,n,indx,r);
        for (int ii=1; ii<=n; ++ii) x(ii) -= r(ii);
    }

    void Inverse(CFARRAY matrix, int nmat)
    {
        using namespace CMLIB;

        matrix.array(nmat,nmat);

        CFARRAY work;
        work.allocate(nmat);
        IFARRAY pivot;
        pivot.allocate(nmat);
        CFARRAY det;
        det.allocate(2);
        int info;
        CGEFA(matrix,nmat,nmat,pivot,info);
        if (info!=0) THROW("Singular matrix in LUdecompose");
        CGEDI(matrix,nmat,nmat,pivot,det,work,1);
    }

    void Inverse(DFARRAY matrix, int nmat)
    {
        using namespace CMLIB;

        matrix.array(nmat,nmat);

        DFARRAY work;
        work.allocate(nmat);
        IFARRAY pivot;
        pivot.allocate(nmat);
        DFARRAY det;
        det.allocate(2);
        int info;
        DGEFA(matrix,nmat,nmat,pivot,info);
        if (info!=0) THROW("Singular matrix in LUdecompose");
        DGEDI(matrix,nmat,nmat,pivot,det,work,1);
    }

    namespace CMLIB {

        static double  CABS1(COMPLEX ZDUM) {
            return fabs(real(ZDUM)) + fabs(imag(ZDUM));
        }

        template <class T> inline T mIn(const T a,const T b) {
            return (a>b)? b : a;
        }
        template <class T> inline T mAx(const T a,const T b) {
            return (a<b)? b : a;
        }

        void CGEEV(DFARRAY A,int LDA, int N, DFARRAY E,
                   DFARRAY V, int LDV, DFARRAY WORK, int JOB, int& INFO)
        {





































































































            int I,IHI,ILO,J,K,L,MDIM,M;
            A.array(2*LDA,N);
            E.array(2*N);
            V.array(2*LDV,N);
            WORK.array(3*N);


            if (N > LDA) THROW("CGEEV-N  >  LDA.");
            if (N < 1 ) THROW("CGEEV-N  <  1");
            if (N == 1 && JOB == 0) goto L35;
            MDIM = 2 * LDA;

            if (JOB == 0) goto L5;
            if (N > LDV) THROW("CGEEV-JOB != 0, AND N  >  LDV.");
            if (N == 1) goto L35;



            MDIM = mIn(MDIM,2 * LDV);
            if (LDA < LDV) THROW("CGEEV-LDA < LDV,  ELEMENTS OF V OTHER THAN THE N BY N OUTPUT ELEMENTS HAVE BEEN CHANGED.");
            if (LDA <= LDV) goto L5;
            SCATMECH_output << "CGEEV-LDA > LDV, ELEMENTS OF A OTHER THAN THE N BY N INPUT ELEMENTS HAVE BEEN CHANGED." << endl;
            L = N - 1;
            for (J=1; J<=L; ++J) {
                I = 2 * N;
                M = 1+J*2*LDV;
                K = 1+J*2*LDA;
                SCOPY(I,A(K),1,A(M),1);
            }
L5:



            for (J=1; J<=N; ++J) {
                K = (J-1) * MDIM +1;
                L = K + N;
                SCOPY(N,A(K+1),2,WORK(1),1);
                SCOPY(N,A(K),2,A(K),1);
                SCOPY(N,WORK(1),1,A(L),1);
            }



            CBAL(MDIM,N,A(1),A(N+1),ILO,IHI,WORK(1));
            CORTH(MDIM,N,ILO,IHI,A(1),A(N+1),WORK(N+1),WORK(2*N+1));
            if (JOB != 0) goto L10;



            COMQR(MDIM,N,ILO,IHI,A(1),A(N+1),E(1),E(N+1),INFO);
            goto L30;



L10:
            COMQR2(MDIM,N,ILO,IHI,WORK(N+1),WORK(2*N+1),A(1),A(N+1),E(1),E(N+1),V(1),V(N+1),INFO);
            if (INFO != 0) goto L30;
            CBABK2(MDIM,N,ILO,IHI,WORK(1),N,V(1),V(N+1));



            for (J=1; J<=N; ++J) {
                K = (J-1) * MDIM + 1;
                I = (J-1) * 2 * LDV + 1;
                L = K + N;
                SCOPY(N,V(K),1,WORK(1),1);
                SCOPY(N,V(L),1,V(I+1),2);
                SCOPY(N,WORK(1),1,V(I),2);
            }



L30:
            SCOPY(N,E(1),1,WORK(1),1);
            SCOPY(N,E(N+1),1,E(2),2);
            SCOPY(N,WORK(1),1,E(1),2);
            return;



L35:
            E(1) = A(1);
            E(2) = A(2);
            INFO = 0;
            if (JOB == 0) return;
            V(1) = A(1);
            V(2) = A(2);
            return;
        }


        void SCOPY(int N,DFARRAY SX,int INCX,DFARRAY SY,int INCY)
        {









































            #ifdef USING_BLAS
            dcopy_(&N,(double*)&SX(1),&INCX,(double*)&SY(1),&INCY);
            #else
            int IX,IY,MP1,NS,I,M;

            if (N == 0) return;
            if (INCX == INCY) {
                if (INCX-1<0) goto L5;
                if (INCX-1==0) goto L20;
                else goto L60;
            }
L5:



            IX = 1;
            IY = 1;
            if (INCX < 0) IX = (-N+1)*INCX + 1;
            if (INCY < 0) IY = (-N+1)*INCY + 1;
            for (I=1; I<=N; ++I) {
                SY[IY-1] = SX[IX-1];
                IX = IX + INCX;
                IY = IY + INCY;
            }
            return;






L20:
            M = N%7;
            if ( M  ==  0 ) goto L40;
            for (I=1; I<=M; ++I) {
                SY[I-1] = SX[I-1];
            }
            if ( N  <  7 ) return;
L40:
            MP1 = M + 1;
            for (I=MP1; I<=N; I+=7) {
                SY[I-1] = SX[I-1];
                SY[I] = SX[I];
                SY[I + 1] = SX[I + 1];
                SY[I + 2] = SX[I + 2];
                SY[I + 3] = SX[I + 3];
                SY[I + 4] = SX[I + 4];
                SY[I + 5] = SX[I + 5];
            }
            return;



L60:
            NS = N*INCX;
            for (I=1; I<=NS; I+=INCX) {
                SY[I-1] = SX[I-1];
            }
            return;
            #endif
        }

        void ZCOPY(int N,CFARRAY ZX,int INCX,CFARRAY ZY,int INCY)
        {










            if (N<0) return;
            if (INCX==1 && INCY==1) {
                for (int I=1; I<=N; ++I) {
                    ZY(I) = ZX(I);
                }
            } else {




                int IX = 1;
                int IY = 1;
                if (INCX<0) IX = (-N+1)*INCX + 1;
                if (INCY<0) IY = (-N+1)*INCY + 1;
                for (int I=1; I<=N; ++I) {
                    ZY(IY) = ZX(IX);
                    IX = IX + INCX;
                    IY = IY + INCY;
                }
            }
        }

        void CBABK2(int& NM,int& N,int& LOW,int& IGH,DFARRAY SCALE,int& M,DFARRAY ZR,DFARRAY ZI)
        {























































            int I,J,K,II;
            ZR.array(NM,M);
            ZI.array(NM,M);
            SCALE.array(N);

            double S;


            if (M  ==  0) goto L200;
            if (IGH  ==  LOW) goto L120;

            for (I=LOW; I<=IGH; ++I) {
                S = SCALE(I);



                for (J=1; J<=M; ++J) {
                    ZR(I,J) = ZR(I,J) * S;
                    ZI(I,J) = ZI(I,J) * S;
                }
            }


L120:
            for (II=1; II<=N; ++II) {
                I = II;
                if (I  >=  LOW  && I  <=  IGH) goto L140;
                if (I  <  LOW) I = LOW - II;
                K = (int)SCALE(I);
                if (K  ==  I) goto L140;

                for (J=1; J<=M; ++J) {
                    S = ZR(I,J);
                    ZR(I,J) = ZR(K,J);
                    ZR(K,J) = S;
                    S = ZI(I,J);
                    ZI(I,J) = ZI(K,J);
                    ZI(K,J) = S;
                }
L140:
                continue;
            }

L200:
            return;
        }

        void CBAL(int& NM,int& N,DFARRAY AR,DFARRAY AI,int& LOW,int& IGH,DFARRAY SCALE)
        {





































































            int I,J,K,L,M,JJ,IEXC;
            AR.array(NM,N);
            AI.array(NM,N);
            SCALE.array(N);
            double C,F,G,R,S,B2,RADIX;
            bool NOCONV;





            RADIX = 16;

            B2 = RADIX * RADIX;
            K = 1;
            L = N;
            goto L100;


L20:
            SCALE(M) = J;
            if (J  ==  M) goto L50;

            for (I=1; I<=L; ++I) {
                F = AR(I,J);
                AR(I,J) = AR(I,M);
                AR(I,M) = F;
                F = AI(I,J);
                AI(I,J) = AI(I,M);
                AI(I,M) = F;
            }

            for (I=K; I<=N; ++I) {
                F = AR(J,I);
                AR(J,I) = AR(M,I);
                AR(M,I) = F;
                F = AI(J,I);
                AI(J,I) = AI(M,I);
                AI(M,I) = F;
            }

L50:
            if (IEXC==1) goto L80;
            if (IEXC==2) goto L130;



L80:
            if (L  ==  1) goto L280;
            L = L - 1;

L100:
            for (JJ=1; JJ<=L; ++JJ) {
                J = L + 1 - JJ;

                for (I=1; I<=L; ++I) {
                    if (I  ==  J) goto L110;
                    if (AR(J,I) != 0.0E0  ||  AI(J,I) != 0.0E0) goto L120;
L110:
                    continue;
                }

                M = L;
                IEXC = 1;
                goto L20;
L120:
                continue;
            }

            goto L140;


L130:
            K = K + 1;

L140:
            for (J=K; J<=L; ++J) {

                for (I=K; I<=L; ++I) {
                    if (I  ==  J) goto L150;
                    if (AR(I,J) != 0.0E0  ||  AI(I,J) != 0.0E0) goto L170;
L150:
                    continue;
                }

                M = K;
                IEXC = 2;
                goto L20;
L170:
                continue;
            }

            for (I=K; I<=L; ++I) {
                SCALE(I) = 1.0E0;
            }

L190:
            NOCONV = false;

            for (I=K; I<=L; ++I) {
                C = 0.0E0;
                R = 0.0E0;

                for (J=K; J<=L; ++J) {
                    if (J  ==  I) goto L200;
                    C = C + fabs(AR(J,I)) + fabs(AI(J,I));
                    R = R + fabs(AR(I,J)) + fabs(AI(I,J));
L200:
                    continue;
                }

                if (C  ==  0.0E0  ||  R  ==  0.0E0) goto L270;
                G = R / RADIX;
                F = 1.0E0;
                S = C + R;
L210:
                if (C  >=  G) goto L220;
                F = F * RADIX;
                C = C * B2;
                goto L210;
L220:
                G = R * RADIX;
L230:
                if (C  <  G) goto L240;
                F = F / RADIX;
                C = C / B2;
                goto L230;

L240:
                if ((C + R) / F  >=  0.95E0 * S) goto L270;
                G = 1.0E0 / F;
                SCALE(I) = SCALE(I) * F;
                NOCONV = true;

                for (J=K; J<=N; ++J) {
                    AR(I,J) = AR(I,J) * G;
                    AI(I,J) = AI(I,J) * G;
                }

                for (J=1; J<=L; ++J) {
                    AR(J,I) = AR(J,I) * F;
                    AI(J,I) = AI(J,I) * F;
                    continue;
                }

L270:
                continue;
            }

            if (NOCONV) goto L190;

L280:
            LOW = K;
            IGH = L;
            return;
        }

        void COMQR(int& NM,int& N,int& LOW,int& IGH,DFARRAY HR,DFARRAY HI,DFARRAY WR,DFARRAY WI,int& IERR)
        {







































































            int I,J,L,EN,LL,ITN,ITS,LP1,ENM1;
            HR.array(NM,N);
            HI.array(NM,N);
            WR.array(N);
            WI.array(N);
            double SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,S1,S2;

            IERR = 0;
            if (LOW  ==  IGH) goto L180;

            L = LOW + 1;

            for (I=L; I<=IGH; ++I) {
                LL = mIn(I+1,IGH);
                if (HI(I,I-1)  ==  0.0E0) goto L170;
                NORM = PYTHAG(HR(I,I-1),HI(I,I-1));
                YR = HR(I,I-1) / NORM;
                YI = HI(I,I-1) / NORM;
                HR(I,I-1) = NORM;
                HI(I,I-1) = 0.0E0;

                for (J=I; J<=IGH; ++J) {
                    SI = YR * HI(I,J) - YI * HR(I,J);
                    HR(I,J) = YR * HR(I,J) + YI * HI(I,J);
                    HI(I,J) = SI;
                    continue;
                }

                for (J=LOW; J<=LL; ++J) {
                    SI = YR * HI(J,I) + YI * HR(J,I);
                    HR(J,I) = YR * HR(J,I) - YI * HI(J,I);
                    HI(J,I) = SI;
                    continue;
                }
L170:
                continue;
            }

L180:
            for (I=1; I<=N; ++I) {
                if (I  >=  LOW  && I  <=  IGH) goto L200;
                WR(I) = HR(I,I);
                WI(I) = HI(I,I);
L200:
                continue;
            }

            EN = IGH;
            TR = 0.0E0;
            TI = 0.0E0;
            ITN = 30*N;

L220:
            if (EN  <  LOW) goto L1001;
            ITS = 0;
            ENM1 = EN - 1;


L240:
            for (LL = LOW; LL<=EN; ++LL) {
                L = EN + LOW - LL;
                if (L  ==  LOW) goto L300;
                S1 = fabs(HR(L-1,L-1)) + fabs(HI(L-1,L-1))+ fabs(HR(L,L)) +fabs(HI(L,L));
                S2 = S1 + fabs(HR(L,L-1));
                if (S2  ==  S1) goto L300;
                continue;
            }

L300:
            if (L  ==  EN) goto L660;
            if (ITN  ==  0) goto L1000;
            if (ITS  ==  10  ||  ITS  ==  20) goto L320;
            SR = HR(EN,EN);
            SI = HI(EN,EN);
            XR = HR(ENM1,EN) * HR(EN,ENM1);
            XI = HI(ENM1,EN) * HR(EN,ENM1);
            if (XR  ==  0.0E0  && XI  ==  0.0E0) goto L340;
            YR = (HR(ENM1,ENM1) - SR) / 2.0E0;
            YI = (HI(ENM1,ENM1) - SI) / 2.0E0;
            CSROOT(sqr(YR)-sqr(YI)+XR,2.0E0*YR*YI+XI,ZZR,ZZI);
            if (YR * ZZR + YI * ZZI  >=  0.0E0) goto L310;
            ZZR = -ZZR;
            ZZI = -ZZI;
L310:
            CDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI);
            SR = SR - XR;
            SI = SI - XI;
            goto L340;

L320:
            SR = fabs(HR(EN,ENM1)) + fabs(HR(ENM1,EN-2));
            SI = 0.0E0;

L340:
            for (I=LOW; I<=EN; ++I) {
                HR(I,I) = HR(I,I) - SR;
                HI(I,I) = HI(I,I) - SI;
                continue;
            }

            TR = TR + SR;
            TI = TI + SI;
            ITS = ITS + 1;
            ITN = ITN - 1;

            LP1 = L + 1;

            for (I=LP1; I<=EN; ++I) {
                SR = HR(I,I-1);
                HR(I,I-1) = 0.0E0;
                NORM = PYTHAG(PYTHAG(HR(I-1,I-1),HI(I-1,I-1)),SR);
                XR = HR(I-1,I-1) / NORM;
                WR(I-1) = XR;
                XI = HI(I-1,I-1) / NORM;
                WI(I-1) = XI;
                HR(I-1,I-1) = NORM;
                HI(I-1,I-1) = 0.0E0;
                HI(I,I-1) = SR / NORM;

                for (J=I; J<=EN; ++J) {
                    YR = HR(I-1,J);
                    YI = HI(I-1,J);
                    ZZR = HR(I,J);
                    ZZI = HI(I,J);
                    HR(I-1,J) = XR * YR + XI * YI + HI(I,I-1) * ZZR;
                    HI(I-1,J) = XR * YI - XI * YR + HI(I,I-1) * ZZI;
                    HR(I,J) = XR * ZZR - XI * ZZI - HI(I,I-1) * YR;
                    HI(I,J) = XR * ZZI + XI * ZZR - HI(I,I-1) * YI;
                    continue;
                }
                continue;
            }

            SI = HI(EN,EN);
            if (SI  ==  0.0E0) goto L540;
            NORM = PYTHAG(HR(EN,EN),SI);
            SR = HR(EN,EN) / NORM;
            SI = SI / NORM;
            HR(EN,EN) = NORM;
            HI(EN,EN) = 0.0E0;

L540:
            for (J=LP1; J<=EN; ++J) {
                XR = WR(J-1);
                XI = WI(J-1);

                for (I=L; I<=J; ++I) {
                    YR = HR(I,J-1);
                    YI = 0.0E0;
                    ZZR = HR(I,J);
                    ZZI = HI(I,J);
                    if (I  ==  J) goto L560;
                    YI = HI(I,J-1);
                    HI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI;
L560:
                    HR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR;
                    HR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR;
                    HI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI;
                    continue;
                }

                continue;
            }

            if (SI  ==  0.0E0) goto L240;

            for (I=L; I<=EN; ++I) {
                YR = HR(I,EN);
                YI = HI(I,EN);
                HR(I,EN) = SR * YR - SI * YI;
                HI(I,EN) = SR * YI + SI * YR;
                continue;
            }

            goto L240;

L660:
            WR(EN) = HR(EN,EN) + TR;
            WI(EN) = HI(EN,EN) + TI;
            EN = ENM1;
            goto L220;


L1000:
            IERR = EN;
L1001:
            return;
        }


        void COMQR2(int& NM,int& N,int& LOW,int& IGH,DFARRAY ORTR,DFARRAY ORTI,
                    DFARRAY HR,DFARRAY HI,DFARRAY WR,DFARRAY WI,DFARRAY ZR,DFARRAY ZI,int& IERR)
        {




















































































            int I,J,K,L,M,EN,II,JJ,LL,NN,IP1;
            int ITN,ITS,LP1,ENM1,IEND;
            HR.array(NM,N);
            HI.array(NM,N);
            WR.array(N);
            WI.array(N);
            ZR.array(NM,N);
            ZI.array(NM,N);
            ORTR.array(IGH);
            ORTI.array(IGH);
            double SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,S1,S2;

            IERR = 0;

            for (I=1; I<=N; ++I) {
                for (J=1; J<=N; ++J) {
                    ZR(I,J) = 0.0E0;
                    ZI(I,J) = 0.0E0;
                    if (I  ==  J) ZR(I,J) = 1.0E0;
                    continue;
                }
            }


            IEND = IGH - LOW - 1;
            if (IEND<0) goto L180;
            if (IEND==0) goto L150;
            if (IEND>0) goto L105;

L105:
            for (II=1; II<=IEND; ++II) {
                I = IGH - II;
                if (ORTR(I)  ==  0.0E0  && ORTI(I)  ==  0.0E0) goto L140;
                if (HR(I,I-1)  ==  0.0E0  && HI(I,I-1)  ==  0.0E0) goto L140;

                NORM = HR(I,I-1) * ORTR(I) + HI(I,I-1) * ORTI(I);
                IP1 = I + 1;

                for (K=IP1; K<=IGH; ++K) {
                    ORTR(K) = HR(K,I-1);
                    ORTI(K) = HI(K,I-1);
                    continue;
                }

                for (J=I; J<=IGH; ++J) {
                    SR = 0.0E0;
                    SI = 0.0E0;

                    for (K=I; K<=IGH; ++K) {
                        SR = SR + ORTR(K) * ZR(K,J) + ORTI(K) * ZI(K,J);
                        SI = SI + ORTR(K) * ZI(K,J) - ORTI(K) * ZR(K,J);
                        continue;
                    }

                    SR = SR / NORM;
                    SI = SI / NORM;

                    for (K=I; K<=IGH; ++K) {
                        ZR(K,J) = ZR(K,J) + SR * ORTR(K) - SI * ORTI(K);
                        ZI(K,J) = ZI(K,J) + SR * ORTI(K) + SI * ORTR(K);
                        continue;
                    }

                    continue;
                }

L140:
                continue;
            }

L150:
            L = LOW + 1;

            for (I=L; I<=IGH; ++I) {
                LL = mIn(I+1,IGH);
                if (HI(I,I-1)  ==  0.0E0) goto L170;
                NORM = PYTHAG(HR(I,I-1),HI(I,I-1));
                YR = HR(I,I-1) / NORM;
                YI = HI(I,I-1) / NORM;
                HR(I,I-1) = NORM;
                HI(I,I-1) = 0.0E0;

                for (J=I; J<=N; ++J) {
                    SI = YR * HI(I,J) - YI * HR(I,J);
                    HR(I,J) = YR * HR(I,J) + YI * HI(I,J);
                    HI(I,J) = SI;
                    continue;
                }

                for (J=1; J<=LL; ++J) {
                    SI = YR * HI(J,I) + YI * HR(J,I);
                    HR(J,I) = YR * HR(J,I) - YI * HI(J,I);
                    HI(J,I) = SI;
                    continue;
                }

                for (J=LOW; J<=IGH; ++J) {
                    SI = YR * ZI(J,I) + YI * ZR(J,I);
                    ZR(J,I) = YR * ZR(J,I) - YI * ZI(J,I);
                    ZI(J,I) = SI;
                    continue;
                }

L170:
                continue;
            }

L180:
            for (I=1; I<=N; ++I) {
                if (I  >=  LOW  && I  <=  IGH) goto L200;
                WR(I) = HR(I,I);
                WI(I) = HI(I,I);
L200:
                continue;
            }

            EN = IGH;
            TR = 0.0E0;
            TI = 0.0E0;
            ITN = 30*N;

L220:
            if (EN  <  LOW) goto L680;
            ITS = 0;
            ENM1 = EN - 1;


L240:
            for (LL=LOW; LL<=EN; ++LL) {
                L = EN + LOW - LL;
                if (L  ==  LOW) goto L300;
                S1 = fabs(HR(L-1,L-1)) + fabs(HI(L-1,L-1)) + fabs(HR(L,L)) +fabs(HI(L,L));
                S2 = S1 + fabs(HR(L,L-1));
                if (S2  ==  S1) goto L300;
                continue;
            }

L300:
            if (L  ==  EN) goto L660;
            if (ITN  ==  0) goto L1000;
            if (ITS  ==  10  ||  ITS  ==  20) goto L320;
            SR = HR(EN,EN);
            SI = HI(EN,EN);
            XR = HR(ENM1,EN) * HR(EN,ENM1);
            XI = HI(ENM1,EN) * HR(EN,ENM1);
            if (XR  ==  0.0E0  && XI  ==  0.0E0) goto L340;
            YR = (HR(ENM1,ENM1) - SR) / 2.0E0;
            YI = (HI(ENM1,ENM1) - SI) / 2.0E0;
            CSROOT(sqr(YR)-sqr(YI)+XR,2.0E0*YR*YI+XI,ZZR,ZZI);
            if (YR * ZZR + YI * ZZI  >=  0.0E0) goto L310;
            ZZR = -ZZR;
            ZZI = -ZZI;
L310:
            CDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI);
            SR = SR - XR;
            SI = SI - XI;
            goto L340;

L320:
            SR = fabs(HR(EN,ENM1)) + fabs(HR(ENM1,EN-2));
            SI = 0.0E0;

L340:
            for (I=LOW; I<=EN; ++I) {
                HR(I,I) = HR(I,I) - SR;
                HI(I,I) = HI(I,I) - SI;
                continue;
            }

            TR = TR + SR;
            TI = TI + SI;
            ITS = ITS + 1;
            ITN = ITN - 1;

            LP1 = L + 1;

            for (I=LP1; I<=EN; ++I) {
                SR = HR(I,I-1);
                HR(I,I-1) = 0.0E0;
                NORM = PYTHAG(PYTHAG(HR(I-1,I-1),HI(I-1,I-1)),SR);
                XR = HR(I-1,I-1) / NORM;
                WR(I-1) = XR;
                XI = HI(I-1,I-1) / NORM;
                WI(I-1) = XI;
                HR(I-1,I-1) = NORM;
                HI(I-1,I-1) = 0.0E0;
                HI(I,I-1) = SR / NORM;

                for (J=I; J<=N; ++J) {
                    YR = HR(I-1,J);
                    YI = HI(I-1,J);
                    ZZR = HR(I,J);
                    ZZI = HI(I,J);
                    HR(I-1,J) = XR * YR + XI * YI + HI(I,I-1) * ZZR;
                    HI(I-1,J) = XR * YI - XI * YR + HI(I,I-1) * ZZI;
                    HR(I,J) = XR * ZZR - XI * ZZI - HI(I,I-1) * YR;
                    HI(I,J) = XR * ZZI + XI * ZZR - HI(I,I-1) * YI;
                    continue;
                }
                continue;
            }

            SI = HI(EN,EN);
            if (SI  ==  0.0E0) goto L540;
            NORM = PYTHAG(HR(EN,EN),SI);
            SR = HR(EN,EN) / NORM;
            SI = SI / NORM;
            HR(EN,EN) = NORM;
            HI(EN,EN) = 0.0E0;
            if (EN  ==  N) goto L540;
            IP1 = EN + 1;

            for (J=IP1; J<=N; ++J) {
                YR = HR(EN,J);
                YI = HI(EN,J);
                HR(EN,J) = SR * YR + SI * YI;
                HI(EN,J) = SR * YI - SI * YR;
                continue;
            }


L540:
            for (J=LP1; J<=EN; ++J) {
                XR = WR(J-1);
                XI = WI(J-1);

                for (I=1; I<=J; ++I) {
                    YR = HR(I,J-1);
                    YI = 0.0E0;
                    ZZR = HR(I,J);
                    ZZI = HI(I,J);
                    if (I  ==  J) goto L560;
                    YI = HI(I,J-1);
                    HI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI;
L560:
                    HR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR;
                    HR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR;
                    HI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI;
                    continue;
                }

                for (I=LOW; I<=IGH; ++I) {
                    YR = ZR(I,J-1);
                    YI = ZI(I,J-1);
                    ZZR = ZR(I,J);
                    ZZI = ZI(I,J);
                    ZR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR;
                    ZI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI;
                    ZR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR;
                    ZI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI;
                    continue;
                }

                continue;
            }

            if (SI  ==  0.0E0) goto L240;

            for (I=1; I<=EN; ++I) {
                YR = HR(I,EN);
                YI = HI(I,EN);
                HR(I,EN) = SR * YR - SI * YI;
                HI(I,EN) = SR * YI + SI * YR;
                continue;
            }

            for (I=LOW; I<=IGH; ++I) {
                YR = ZR(I,EN);
                YI = ZI(I,EN);
                ZR(I,EN) = SR * YR - SI * YI;
                ZI(I,EN) = SR * YI + SI * YR;
                continue;
            }

            goto L240;

L660:
            HR(EN,EN) = HR(EN,EN) + TR;
            WR(EN) = HR(EN,EN);
            HI(EN,EN) = HI(EN,EN) + TI;
            WI(EN) = HI(EN,EN);
            EN = ENM1;
            goto L220;


L680:
            NORM = 0.0E0;

            for (I=1; I<=N; ++I) {
                for (J=I; J<=N; ++J) {
                    NORM = NORM + fabs(HR(I,J)) + fabs(HI(I,J));
                    continue;
                }
            }

            if (N  ==  1  ||  NORM  ==  0.0E0) goto L1001;

            for (NN=2; NN<=N; ++NN) {
                EN = N + 2 - NN;
                XR = WR(EN);
                XI = WI(EN);
                ENM1 = EN - 1;

                for (II=1; II<=ENM1; ++II) {
                    I = EN - II;
                    ZZR = HR(I,EN);
                    ZZI = HI(I,EN);
                    if (I  ==  ENM1) goto L760;
                    IP1 = I + 1;

                    for (J=IP1; J<=ENM1; ++J) {
                        ZZR = ZZR + HR(I,J) * HR(J,EN) - HI(I,J) * HI(J,EN);
                        ZZI = ZZI + HR(I,J) * HI(J,EN) + HI(I,J) * HR(J,EN);
                        continue;
                    }

L760:
                    YR = XR - WR(I);
                    YI = XI - WI(I);
                    if (YR != 0.0E0  ||  YI != 0.0E0) goto L775;
                    YR = NORM;
L770:
                    YR = 0.5E0*YR;
                    if (NORM + YR  >  NORM) goto L770;
                    YR = 2.0E0*YR;
L775:
                    CDIV(ZZR,ZZI,YR,YI,HR(I,EN),HI(I,EN));
                    continue;
                }
                continue;
            }

            ENM1 = N - 1;

            for (I=1; I<=ENM1; ++I) {
                if (I  >=  LOW  && I  <=  IGH) goto L840;
                IP1 = I + 1;

                for (J=IP1; J<=N; ++J) {
                    ZR(I,J) = HR(I,J);
                    ZI(I,J) = HI(I,J);
                    continue;
                }
L840:
                continue;
            }



            for (JJ=LOW; JJ<=ENM1; ++JJ) {
                J = N + LOW - JJ;
                M = mIn(J-1,IGH);

                for (I=LOW; I<=IGH; ++I) {
                    ZZR = ZR(I,J);
                    ZZI = ZI(I,J);

                    for (K=LOW; K<=M; ++K) {
                        ZZR = ZZR + ZR(I,K) * HR(K,J) - ZI(I,K) * HI(K,J);
                        ZZI = ZZI + ZR(I,K) * HI(K,J) + ZI(I,K) * HR(K,J);
                        continue;
                    }

                    ZR(I,J) = ZZR;
                    ZI(I,J) = ZZI;
                    continue;
                }
            }
            goto L1001;



L1000:
            IERR = EN;
L1001:
            return;
        }


        void CORTH(int& NM,int& N,int& LOW,int& IGH,DFARRAY AR,DFARRAY AI,DFARRAY ORTR,DFARRAY ORTI)
        {



























































            int I,J,M,II,JJ,LA,MP,KP1;
            AR.array(NM,N);
            AI.array(NM,N);
            ORTR.array(IGH);
            ORTI.array(IGH);
            double F,G,H,FI,FR,SCALE;



            LA = IGH - 1;
            KP1 = LOW + 1;
            if (LA  <  KP1) goto L200;

            for (M=KP1; M<=LA; ++M) {
                H = 0.0E0;
                ORTR(M) = 0.0E0;
                ORTI(M) = 0.0E0;
                SCALE = 0.0E0;

                for (I=M; I<=IGH; ++I) {
                    SCALE = SCALE + fabs(AR(I,M-1)) + fabs(AI(I,M-1));
                }

                if (SCALE  ==  0.0E0) goto L180;
                MP = M + IGH;

                for (II=M; II<=IGH; ++II) {
                    I = MP - II;
                    ORTR(I) = AR(I,M-1) / SCALE;
                    ORTI(I) = AI(I,M-1) / SCALE;
                    H = H + ORTR(I) * ORTR(I) + ORTI(I) * ORTI(I);
                    continue;
                }

                G = sqrt(H);
                F = PYTHAG(ORTR(M),ORTI(M));
                if (F  ==  0.0E0) goto L103;
                H = H + F * G;
                G = G / F;
                ORTR(M) = (1.0E0 + G) * ORTR(M);
                ORTI(M) = (1.0E0 + G) * ORTI(M);
                goto L105;

L103:
                ORTR(M) = G;
                AR(M,M-1) = SCALE;

L105:
                for (J=M; J<=N; ++J) {
                    FR = 0.0E0;
                    FI = 0.0E0;

                    for (II=M; II<=IGH; ++II) {
                        I = MP - II;
                        FR = FR + ORTR(I) * AR(I,J) + ORTI(I) * AI(I,J);
                        FI = FI + ORTR(I) * AI(I,J) - ORTI(I) * AR(I,J);
                        continue;
                    }

                    FR = FR / H;
                    FI = FI / H;

                    for (I=M; I<=IGH; ++I) {
                        AR(I,J) = AR(I,J) - FR * ORTR(I) + FI * ORTI(I);
                        AI(I,J) = AI(I,J) - FR * ORTI(I) - FI * ORTR(I);
                        continue;
                    }

                    continue;
                }

                for (I=1; I<=IGH; ++I) {
                    FR = 0.0E0;
                    FI = 0.0E0;

                    for (JJ=M; JJ<=IGH; ++JJ) {
                        J = MP - JJ;
                        FR = FR + ORTR(J) * AR(I,J) - ORTI(J) * AI(I,J);
                        FI = FI + ORTR(J) * AI(I,J) + ORTI(J) * AR(I,J);
                        continue;
                    }

                    FR = FR / H;
                    FI = FI / H;

                    for (J=M; J<=IGH; ++J) {
                        AR(I,J) = AR(I,J) - FR * ORTR(J) - FI * ORTI(J);
                        AI(I,J) = AI(I,J) + FR * ORTI(J) - FI * ORTR(J);
                        continue;
                    }

                    continue;
                }

                ORTR(M) = SCALE * ORTR(M);
                ORTI(M) = SCALE * ORTI(M);
                AR(M,M-1) = -G * AR(M,M-1);
                AI(M,M-1) = -G * AI(M,M-1);
L180:
                continue;
            }

L200:
            return;
        }


        void CSROOT(const double& XR, const double& XI, DFARRAY YR, DFARRAY YI)
        {







            double S,TR,TI;



            TR = XR;
            TI = XI;
            S = sqrt(0.5E0*(PYTHAG(TR,TI) + fabs(TR)));
            if (TR  >=  0.0E0) YR = S;
            if (TI  <  0.0E0) S = -S;
            if (TR  <=  0.0E0) YI = S;
            if (TR  <  0.0E0) YR = 0.5E0*(TI/YI);
            if (TR  >  0.0E0) YI = 0.5E0*(TI/YR);
            return;
        }

        double PYTHAG(const double& A,const double& B)
        {







            double P,Q,R,S,T;

            P = mAx(fabs(A),fabs(B));
            Q = mIn(fabs(A),fabs(B));
            if (Q  ==  0.0E0) goto L20;
L10:
            R = sqr(Q/P);
            T = 4.0E0 + R;
            if (T  ==  4.0E0) goto L20;
            S = R/T;
            P = P + 2.0E0*P*S;
            Q = Q*S;
            goto L10;
L20:
            return (P);
        }


        void CDIV(const double& AR,const double& AI,const double& BR,const double& BI,DFARRAY CR,DFARRAY CI)
        {







            double S,ARS,AIS,BRS,BIS;

            S = fabs(BR) + fabs(BI);
            ARS = AR/S;
            AIS = AI/S;
            BRS = BR/S;
            BIS = BI/S;
            S = sqr(BRS) + sqr(BIS);
            CR = (ARS*BRS + AIS*BIS)/S;
            CI = (AIS*BRS - ARS*BIS)/S;
            return;
        }

        void CGEFA(FARRAY<COMPLEX> A, int LDA, int N, FARRAY<int> IPVT, int& INFO)
        {



























































            A.array(LDA,N);
            IPVT.array(N);

            COMPLEX T;
            int J,K,KP1,L,NM1;
            COMPLEX ZDUM;

            INFO = 0;
            NM1 = N - 1;
            if (NM1 >= 1) {
                for (K = 1; K<=NM1; ++K) {
                    KP1 = K + 1;
                    L = ICAMAX(N-K+1,A(K,K),1) + K - 1;

                    IPVT(K) = L;
                    if (CABS1(A(L,K)) == 0.0E0) {
                        INFO = K;
                    } else {
                        if (L != K) {
                            T = A(L,K);
                            A(L,K) = A(K,K);
                            A(K,K) = T;
                        }
                        T = -COMPLEX(1.0E0,0.0E0)/A(K,K);
                        CSCAL(N-K,T,A(K+1,K),1);
                        for (J=KP1; J<=N; ++J) {
                            T = A(L,J);
                            if (L != K) {
                                A(L,J) = A(K,J);
                                A(K,J) = T;
                            }
                            CAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1);
                        }
                    }
                }
            }
            IPVT(N) = N;
            if (CABS1(A(N,N)) == 0.0E0) INFO = N;
            return;
        }

        void CAXPY(int N, COMPLEX CA, FARRAY<COMPLEX> CX, int INCX, FARRAY<COMPLEX> CY, int INCY)
        {










































            #ifdef USING_BLAS
            zaxpy_(&N,(cmplx*)&CA,(cmplx*)&CX(1),&INCX,(cmplx*)&CY(1),&INCY);
            #else
            CX.array(N);
            CY.array(N);
            int KX,KY,I;
            COMPLEX CANORM = fabs(real(CA)) + fabs(imag(CA));
            if (N<=0||CANORM==0.E0) return;
            if (INCX==INCY&&INCX>0) {
                int NS = N*INCX;
                for (I=1; I<=NS; I+=INCX) {
                    CY(I) = CA*CX(I) + CY(I);
                }
                return;
            }
            KX = 1;
            KY = 1;
            if (INCX<0) KX = 1+(1-N)*INCX;
            if (INCY<0) KY = 1+(1-N)*INCY;
            for (I=1; I<=N; ++I) {
                CY(KY) = CY(KY) + CA*CX(KX);
                KX = KX + INCX;
                KY = KY + INCY;
            }
            return;
            #endif
        }


        void CSCAL(int N, COMPLEX CA, FARRAY<COMPLEX> CX, int INCX)
        {






































            #ifdef USING_BLAS
            zscal_(&N, (cmplx*)&CA, (cmplx*)&CX(1), &INCX);
            #else
            CX.array(N);
            if (N <= 0) return;
            int NS = N*INCX;
            for (int I=1; I<=NS; I+=INCX) {
                CX(I) = CA*CX(I);
            }
            return;
            #endif
        }

        int ICAMAX(int N, FARRAY<COMPLEX> CX, int INCX)
        {







































            #ifdef USING_BLAS
            return izamax_(&N,(cmplx*)&CX(1),&INCX);
            #else
            CX.array(N);
            int _ICAMAX = 0;
            if (N<=0) return _ICAMAX;
            _ICAMAX = 1;
            if (N <= 1) return _ICAMAX;
            int NS = N*INCX;
            int II = 1;
            double SUMMAX = fabs(real(CX(1))) + fabs(imag(CX(1)));
            for (int I=1; I<=NS; I+=INCX) {
                double SUMRI = fabs(real(CX(I))) + fabs(imag(CX(I)));
                if (SUMMAX<SUMRI) {
                    SUMMAX = SUMRI;
                    _ICAMAX = II;
                }
                II = II + 1;
            }
            return _ICAMAX;
            #endif
        }

        void CGESL(FARRAY<COMPLEX> A, int LDA, int N, FARRAY<int> IPVT, FARRAY<COMPLEX> B, int JOB)
        {







































































            A.array(LDA,N);
            B.array(N);
            IPVT.array(N);

            COMPLEX T;
            int K,KB,L,NM1;
            NM1 = N - 1;
            if (JOB == 0) {




                if (NM1 >= 1) {
                    for (K=1; K<=NM1; ++K) {
                        L = IPVT(K);
                        T = B(L);
                        if (L != K) {
                            B(L) = B(K);
                            B(K) = T;
                        }
                        CAXPY(N-K,T,A(K+1,K),1,B(K+1),1);
                    }




                    for (KB=1; KB<=N; ++KB) {
                        K = N + 1 - KB;
                        B(K) = B(K)/A(K,K);
                        T = -B(K);
                        CAXPY(K-1,T,A(1,K),1,B(1),1);
                    }
                }
            }   else {





                for (K=1; K<=N; ++K) {
                    T = CDOTC(K-1,A(1,K),1,B(1),1);
                    B(K) = (B(K) - T)/conj(A(K,K));
                }



                if (NM1 >= 1) {
                    for (KB=1; KB<=NM1; ++KB) {
                        K = N - KB;
                        B(K) = B(K) + CDOTC(N-K,A(K+1,K),1,B(K+1),1);
                        L = IPVT(K);
                        if (L != K) {
                            T = B(L);
                            B(L) = B(K);
                            B(K) = T;
                        }
                    }
                }
            }
        }

        COMPLEX CDOTC(int N, FARRAY<COMPLEX> CX, int INCX, FARRAY<COMPLEX> CY, int INCY)
        {










































            #ifdef USING_BLAS
            cmplx result = zdotc_(&N,(cmplx*)&CX(1),&INCX,(cmplx*)&CY(1),&INCY);
            return COMPLEX(result.re,result.im);
            #else
            COMPLEX _CDOTC(0.,0.);
            int I,KX,KY;
            if (N < 0) return _CDOTC;
            if (INCX==INCY&&INCX>0) {
                int NS = N*INCX;
                for (I=1; I<=NS; I+=INCX) {
                    _CDOTC = conj(CX(I))*CY(I) + _CDOTC;
                }
                return _CDOTC;
            } else {
                KX = 1;
                KY = 1;
                if (INCX<0) KX = 1+(1-N)*INCX;
                if (INCY<0) KY = 1+(1-N)*INCY;
                for (I=1; I<=N; ++I) {
                    _CDOTC = _CDOTC + conj(CX(KX))*CY(KY);
                    KX = KX + INCX;
                    KY = KY + INCY;
                }
                return _CDOTC;
            }
            #endif
        }

        void CGEDI(CFARRAY A, int LDA, int N, IFARRAY IPVT, CFARRAY DET, CFARRAY WORK, int JOB)
        {





































































            A.array(LDA,N);
            WORK.array(N);
            IPVT.array(N);


            COMPLEX T;
            double TEN;
            int I,J,K,KB,KP1,L,NM1;
            COMPLEX ZDUM;




            if (JOB/10 != 0) {
                DET(1) = COMPLEX(1.0E0,0.0E0);
                DET(2) = COMPLEX(0.0E0,0.0E0);
                TEN = 10.0E0;
                for (I=1; I<=N; ++I) {
                    if (IPVT(I) != I) DET(1) = -DET(1);
                    DET(1) = A(I,I)*DET(1);

                    if (CABS1(DET(1)) == 0.0E0) break;
                    while (CABS1(DET(1)) < 1.0E0) {
                        DET(1) = COMPLEX(TEN,0.0E0)*DET(1);
                        DET(2) = DET(2) - COMPLEX(1.0E0,0.0E0);
                    }
                    while (CABS1(DET(1)) >= TEN) {
                        DET(1) = DET(1)/COMPLEX(TEN,0.0E0);
                        DET(2) = DET(2) + COMPLEX(1.0E0,0.0E0);
                    }
                }
            }



            if (JOB%10 != 0) {
                for (K=1; K<=N; ++K) {
                    A(K,K) = COMPLEX(1.0E0,0.0E0)/A(K,K);
                    T = -A(K,K);
                    CSCAL(K-1,T,A(1,K),1);
                    KP1 = K + 1;
                    if (N >= KP1) {
                        for (J = KP1; J<=N; ++J) {
                            T = A(K,J);
                            A(K,J) = COMPLEX(0.0E0,0.0E0);
                            CAXPY(K,T,A(1,K),1,A(1,J),1);
                        }
                    }
                }



                NM1 = N - 1;
                if (NM1 >= 1) {
                    for (KB=1; KB<=NM1; ++KB) {
                        K = N - KB;
                        KP1 = K + 1;
                        for (I=KP1; I<=N; ++I) {
                            WORK(I) = A(I,K);
                            A(I,K) = COMPLEX(0.0E0,0.0E0);
                        }
                        for (J=KP1; J<=N; ++J) {
                            T = WORK(J);
                            CAXPY(N,T,A(1,J),1,A(1,K),1);
                        }
                        L = IPVT(K);
                        if (L != K) CSWAP(N,A(1,K),1,A(1,L),1);
                    }
                }
            }
        }


        void CSWAP(int N,CFARRAY CX, int INCX, CFARRAY CY, int INCY)
        {








































            #ifdef USING_BLAS
            zswap_(&N,(cmplx*)&CX(1),&INCX,(cmplx*)&CY(1),&INCY);
            #else
            CX.array(N);
            CY.array(N);
            COMPLEX CTEMP;
            int KX,KY,I,NS;

            if (N <= 0) return;
            if (INCX==INCY&&INCX>0) {
                NS = N*INCX;
                for (I=1; I<=NS; I+=INCX) {
                    CTEMP = CX(I);
                    CX(I) = CY(I);
                    CY(I) = CTEMP;
                }

            } else {
                KX = 1;
                KY = 1;
                if (INCX<0) KX = 1+(1-N)*INCX;
                if (INCY<0) KY = 1+(1-N)*INCY;
                for (I=1; I<=N; ++I) {
                    CTEMP = CX(KX);
                    CX(KX) = CY(KY);
                    CY(KY) = CTEMP;
                    KX = KX + INCX;
                    KY = KY + INCY;
                }
            }
            #endif
        }

        bool lsame(char a,char b)
        {
            if (tolower(a)==tolower(b)) return true;
            else return false;
        }

        void ZGEMV(const char *TRANS,int M,int N,COMPLEX ALPHA,CFARRAY A,int LDA,CFARRAY X,int INCX,COMPLEX BETA,CFARRAY Y,int INCY)
        {
            #ifdef USING_BLAS
            zgemv_(TRANS,&M,&N,(cmplx*)&ALPHA,(cmplx*)&A(1),&LDA,(cmplx*)&X(1),&INCX,(cmplx*)&BETA,(cmplx*)&Y(1),&INCY);
            #else
            A.array(LDA,N);
































































































            const COMPLEX ONE(1,0);
            const COMPLEX ZERO(0,0);


            COMPLEX TEMP;
            int INFO,IX,IY,JX,JY,KX,KY,LENX,LENY;
            bool NOCONJ;







            INFO = 0;
            if (!lsame(*TRANS,'N')&&!lsame(*TRANS,'T')&&!lsame(*TRANS,'C')) INFO=1;
            else if (M<0) INFO=2;
            else if (N<0) INFO=3;
            else if (LDA<(1>M?1:M)) INFO=6;
            else if (INCX==0) INFO=8;
            else if (INCY==0) INFO=11;
            if (INFO!=0) throw SCATMECH_exception("Invalid arguments to ZGEMV...Error code "+to_string(INFO));





            if (M==0 || N==0 || (ALPHA==ZERO && BETA==ONE)) return;

            NOCONJ = lsame(*TRANS,'T');




            if (lsame(*TRANS,'N')) {
                LENX = N;
                LENY = M;
            } else {
                LENX = M;
                LENY = N;
            }
            if (INCX>0) {
                KX = 1;
            } else {
                KX = 1 - (LENX-1)*INCX;
            }
            if (INCY>0) {
                KY = 1;
            } else {
                KY = 1 - (LENY-1)*INCY;
            }






            if (BETA!=ONE) {
                if (INCY==1) {
                    if (BETA==ZERO) {
                        for (int I=1; I<=LENY; ++I) {
                            Y(I) = ZERO;
                        }
                    } else {
                        for (int I=1; I<=LENY; ++I) {
                            Y(I) = BETA*Y(I);
                        }
                    }
                } else {
                    IY = KY;
                    if (BETA==ZERO) {
                        for (int I=1; I<=LENY; ++I) {
                            Y(IY) = ZERO;
                            IY = IY + INCY;
                        }
                    } else {
                        for (int I=1; I<=LENY; ++I) {
                            Y(IY) = BETA*Y(IY);
                            IY = IY + INCY;
                        }
                    }
                }
            }
            if (ALPHA==ZERO) return;


            if (lsame(*TRANS,'N')) {



                JX = KX;
                if (INCY==1) {
                    for (int J=1; J<=N; ++J) {
                        if (X(JX)!=ZERO) {
                            TEMP = ALPHA*X(JX);
                            for (int I=1; I<=M; ++I) {
                                Y(I) = Y(I) + TEMP*A(I,J);
                            }
                        }
                        JX = JX + INCX;
                    }
                } else {
                    for (int J=1; J<=N; ++J) {
                        if (X(JX)!=ZERO) {
                            TEMP = ALPHA*X(JX);
                            IY = KY;
                            for (int I=1; I<=M; ++I) {
                                Y(IY) = Y(IY) + TEMP*A(I,J);
                                IY = IY + INCY;
                            }
                        }
                        JX = JX + INCX;
                    }
                }
            } else {



                JY = KY;
                if (INCX==1) {
                    for (int J=1; J<=N; ++J) {
                        TEMP = ZERO;
                        if (NOCONJ) {
                            for (int I=1; I<=M; ++I) {
                                TEMP = TEMP + A(I,J)*X(I);
                            }
                        } else {
                            for (int I=1; I<=M; ++I) {
                                TEMP = TEMP + conj(A(I,J))*X(I);
                            }
                        }
                        Y(JY) = Y(JY) + ALPHA*TEMP;
                        JY = JY + INCY;
                    }
                } else {
                    for (int J=1; J<=N; ++J) {
                        TEMP = ZERO;
                        IX = KX;
                        if (NOCONJ) {
                            for (int I=1; I<=M; ++I) {
                                TEMP = TEMP + A(I,J)*X(IX);
                                IX = IX + INCX;
                            }
                        } else {
                            for (int I=1; I<=M; ++I) {
                                TEMP = TEMP + conj(A(I,J))*X(IX);
                                IX = IX + INCX;
                            }
                        }
                        Y(JY) = Y(JY) + ALPHA*TEMP;
                        JY = JY + INCY;
                    }
                }
            }


            return;



            #endif
        }

        void ZGEMM(const char* TRANSA,const char* TRANSB,int M,int N,int K,COMPLEX ALPHA,CFARRAY A,int LDA,CFARRAY B,int LDB,COMPLEX BETA,CFARRAY C,int LDC)
        {
            #ifdef USING_BLAS
            zgemm_(TRANSA,TRANSB,&M,&N,&K,(cmplx*)&ALPHA,(cmplx*)&A(1),&LDA,(cmplx*)&B(1),&LDB,(cmplx*)&BETA,(cmplx*)&C(1),&LDC);
            #else


            A.array(LDA,K);
            B.array(LDB,K);
            C.array(LDC,N);



























































































































            COMPLEX TEMP;
            int INFO,NCOLA,NROWA,NROWB;
            bool CONJA,CONJB,NOTA,NOTB;

            const COMPLEX ONE(1,0);
            const COMPLEX ZERO(0,0);








            NOTA = lsame(*TRANSA,'N');
            NOTB = lsame(*TRANSB,'N');
            CONJA = lsame(*TRANSA,'C');
            CONJB = lsame(*TRANSB,'C');
            if (NOTA) {
                NROWA = M;
                NCOLA = K;
            } else {
                NROWA = K;
                NCOLA = M;
            }
            if (NOTB) {
                NROWB = K;
            } else {
                NROWB = N;
            }



            INFO = 0;
            if (!NOTA && !CONJA && !lsame(*TRANSA,'T')) INFO = 1;
            else if (!NOTB && !CONJB && !lsame(*TRANSB,'T')) INFO = 2;
            else if (M<0) INFO = 3;
            else if (N<0) INFO = 4;
            else if (K<0) INFO = 5;
            else if (LDA<(1>NROWA?1:NROWA)) INFO = 8;
            else if (LDB<(1>NROWB?1:NROWB)) INFO = 10;
            else if (LDC<(1>M?1:M)) INFO = 13;

            if (INFO!=0) throw SCATMECH_exception("Invalid arguments to ZGEMM...Error code "+to_string(INFO));




            if (M==0 || N==0 || ((ALPHA==ZERO || K==0) && BETA==ONE)) return;



            if (ALPHA==ZERO) {
                if (BETA==ZERO) {
                    for (int J=1; J<=N; ++J) {
                        for (int I=1; I<=M; ++I) {
                            C(I,J) = ZERO;
                        }
                    }
                } else {
                    for (int J=1; J<=N; ++J) {
                        for (int I=1; I<=M; ++I) {
                            C(I,J) = BETA*C(I,J);
                        }
                    }
                }
                return;
            }



            if (NOTB) {
                if (NOTA) {



                    for (int J=1; J<=N; ++J) {
                        if (BETA==ZERO) {
                            for (int I=1; I<=M; ++I) {
                                C(I,J) = ZERO;
                            }
                        } else if (BETA!=ONE) {
                            for (int I=1; I<=M; ++I) {
                                C(I,J) = BETA*C(I,J);
                            }
                        }
                        for (int L=1; L<=K; ++L) {
                            if (B(L,J)!=ZERO) {
                                TEMP = ALPHA*B(L,J);
                                for (int I=1; I<=M; ++I) {
                                    C(I,J) = C(I,J) + TEMP*A(I,L);
                                }
                            }
                        }
                    }
                } else if (CONJA) {



                    for (int J=1; J<=N; ++J) {
                        for (int I=1; I<=M; ++I) {
                            TEMP = ZERO;
                            for (int L=1; L<=K; ++L) {
                                TEMP = TEMP + conj(A(L,I))*B(L,J);
                            }
                            if (BETA==ZERO) {
                                C(I,J) = ALPHA*TEMP;
                            } else {
                                C(I,J) = ALPHA*TEMP + BETA*C(I,J);
                            }
                        }
                    }
                } else {



                    for (int J=1; J<=N; ++J) {
                        for (int I=1; I<=M; ++I) {
                            TEMP = ZERO;
                            for (int L=1; L<=K; ++L) {
                                TEMP = TEMP + A(L,I)*B(L,J);
                            }
                            if (BETA==ZERO) {
                                C(I,J) = ALPHA*TEMP;
                            } else {
                                C(I,J) = ALPHA*TEMP + BETA*C(I,J);
                            }
                        }
                    }
                }
            } else if (NOTA) {
                if (CONJB) {



                    for (int J=1; J<=N; ++J) {
                        if (BETA==ZERO) {
                            for (int I=1; I<=M; ++I) {
                                C(I,J) = ZERO;
                            }
                        } else if (BETA!=ONE) {
                            for (int I=1; I<=M; ++I) {
                                C(I,J) = BETA*C(I,J);
                            }

                        }
                        for (int L=1; L<=K; ++L) {
                            if (B(J,L)!=ZERO) {
                                TEMP = ALPHA*conj(B(J,L));
                                for (int I=1; I<=M; ++I) {
                                    C(I,J) = C(I,J) + TEMP*A(I,L);
                                }
                            }
                        }
                    }
                } else {



                    for (int J=1; J<=N; ++J) {
                        if (BETA==ZERO) {
                            for (int I=1; I<=M; ++I) {
                                C(I,J) = ZERO;
                            }
                        } else if (BETA!=ONE) {
                            for (int I=1; I<=M; ++I) {
                                C(I,J) = BETA*C(I,J);
                            }
                        }
                        for (int L=1; L<=K; ++L) {
                            if (B(J,L)!=ZERO) {
                                TEMP = ALPHA*B(J,L);
                                for (int I=1; I<=M; ++I) {
                                    C(I,J) = C(I,J) + TEMP*A(I,L);
                                }
                            }
                        }
                    }
                }
            } else if (CONJA) {
                if (CONJB) {



                    for (int J=1; J<=N; ++J) {
                        for (int I=1; I<=M; ++I) {
                            TEMP = ZERO;
                            for (int L=1; L<=K; ++L) {
                                TEMP = TEMP + conj(A(L,I))*conj(B(J,L));
                            }
                            if (BETA==ZERO) {
                                C(I,J) = ALPHA*TEMP;
                            } else {
                                C(I,J) = ALPHA*TEMP + BETA*C(I,J);
                            }
                        }
                    }
                } else {



                    for (int J=1; J<=N; ++J) {
                        for (int I=1; I<=M; ++I) {
                            TEMP = ZERO;
                            for (int L=1; L<=K; ++L) {
                                TEMP = TEMP + conj(A(L,I))*B(J,L);
                            }
                            if (BETA==ZERO) {
                                C(I,J) = ALPHA*TEMP;
                            } else {
                                C(I,J) = ALPHA*TEMP + BETA*C(I,J);
                            }
                        }
                    }
                }
            } else {
                if (CONJB) {



                    for (int J=1; J<=N; ++J) {
                        for (int I=1; I<=M; ++I) {
                            TEMP = ZERO;
                            for (int L=1; L<=K; ++L) {
                                TEMP = TEMP + A(L,I)*conj(B(J,L));
                            }
                            if (BETA==ZERO) {
                                C(I,J) = ALPHA*TEMP;
                            } else {
                                C(I,J) = ALPHA*TEMP + BETA*C(I,J);
                            }
                        }
                    }
                } else {



                    for (int J=1; J<=N; ++J) {
                        for (int I=1; I<=M; ++I) {
                            TEMP = ZERO;
                            for (int L=1; L<=K; ++L) {
                                TEMP = TEMP + A(L,I)*B(J,L);
                            }
                            if (BETA==ZERO) {
                                C(I,J) = ALPHA*TEMP;
                            } else {
                                C(I,J) = ALPHA*TEMP + BETA*C(I,J);
                            }
                        }
                    }
                }
            }
            return;



            #endif
        }

        void DGEFA(FARRAY<double> A,int LDA,int N,FARRAY<int> IPVT,int& INFO)
        {
            #ifdef USING_BLAS
            dgefa_((double*)&A(1),&LDA,&N,(int*)&IPVT(1),&INFO);
            #else
























































            A.array(LDA,N);
            double T;
            int J,K,KP1,L,NM1;




            INFO = 0;
            NM1 = N - 1;
            if (NM1 >= 1) {
                for (K=1; K<=NM1; ++K) {
                    KP1 = K + 1;



                    L = IDAMAX(N-K+1,A(K,K),1) + K - 1;
                    IPVT(K) = L;



                    if (A(L,K) != 0.) {



                        if (L!=K) {
                            T = A(L,K);
                            A(L,K) = A(K,K);
                            A(K,K) = T;
                        }



                        T = -1.0/A(K,K);
                        DSCAL(N-K,T,A(K+1,K),1);



                        for (J=KP1; J<=N; ++J) {
                            T = A(L,J);
                            if (L!=K) {
                                A(L,J) = A(K,J);
                                A(K,J) = T;
                            }
                            DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1);
                        }
                    } else {
                        INFO = K;
                    }
                }
            }
            IPVT(N) = N;
            if (A(N,N) == 0.0) INFO = N;
            return;
            #endif
        }

        void DAXPY(int N,double DA,FARRAY<double> DX,int INCX,FARRAY<double> DY,int INCY)
        {
            #ifdef USING_BLAS
            daxpy_(&N,&DA,(double*)&DX(1),&INCX,(double*)&DY(1),&INCY);
            #else








































            DX.array(INCX);
            DY.array(INCY);

            if (N<=0||DA==0.0) return;

            if (INCX!=INCY || INCX-1<0) {



                int IX = 1;
                int IY = 1;
                if (INCX<0) IX = (-N+1)*INCX + 1;
                if (INCY<0) IY = (-N+1)*INCY + 1;
                for (int I=1; I<=N; ++I) {
                    DY(IY) = DY(IY) + DA*DX(IX);
                    IX = IX + INCX;
                    IY = IY + INCY;
                }
                return;
            } else if (INCX-1==0) {






                int M = N%4;
                if ( M != 0 ) {
                    for (int I = 1; I<=M; ++I) {
                        DY(I) = DY(I) + DA*DX(I);
                    }
                    if (N<4) return;
                }
                int MP1 = M + 1;
                for (int I = MP1; I<=N; I+=4) {
                    DY(I) = DY(I) + DA*DX(I);
                    DY(I + 1) = DY(I + 1) + DA*DX(I + 1);
                    DY(I + 2) = DY(I + 2) + DA*DX(I + 2);
                    DY(I + 3) = DY(I + 3) + DA*DX(I + 3);
                }
                return;
            } else {



                int NS = N*INCX;
                for (int I=1; I<=NS; I+=INCX) {
                    DY(I) = DA*DX(I) + DY(I);
                }
                return;
            }
            #endif
        }

        void DSCAL(int N,double DA,FARRAY<double> DX,int INCX)
        {
            #ifdef USING_BLAS
            dscal_(&N,&DA,(double*)&DX(1),&INCX);
            #else




































            DX.array(N*INCX);

            if (N<=0) return;
            if (INCX!=1) {



                int NS = N*INCX;
                for (int I=1; I<=NS; I+=INCX) {
                    DX(I) = DA*DX(I);
                }
                return;
            }






            int M = N%5;
            if (M!=0) {
                for (int I=1; I<=M; ++I) {
                    DX(I) = DA*DX(I);
                }
                if (N<5) return;
            }
            int MP1 = M + 1;
            for (int I=MP1; I<=N; I+=5) {
                DX(I) = DA*DX(I);
                DX(I + 1) = DA*DX(I + 1);
                DX(I + 2) = DA*DX(I + 2);
                DX(I + 3) = DA*DX(I + 3);
                DX(I + 4) = DA*DX(I + 4);
            }
            return;
            #endif
        }

        int IDAMAX(int N,FARRAY<double> DX,int INCX)
        {
            #ifdef USING_BLAS
            idamax_(&N,(double*)&DX(1),&INCX);
            #else




































            double DMAX,XMAG;
            DX.array(N);

            int _IDAMAX = 0;
            if (N<=0) return _IDAMAX;
            _IDAMAX = 1;
            if (N<=1) return _IDAMAX;
            if (INCX!=1) {



                DMAX = fabs(DX(1));
                int NS = N*INCX;
                int II = 1;
                for (int I=1; I<=NS; I+=INCX) {
                    XMAG = fabs(DX(I));
                    if (XMAG>DMAX) {
                        _IDAMAX = II;
                        DMAX = XMAG;
                    }
                    II = II + 1;
                }
                return _IDAMAX;
            }



            DMAX = fabs(DX(1));
            for (int I=2; I<=N; ++I) {
                XMAG = fabs(DX(I));
                if (XMAG>DMAX) {
                    _IDAMAX = I;
                    DMAX = XMAG;
                }
            }
            return _IDAMAX;
            #endif
        }

        void DGEDI(FARRAY<double> A,int LDA,int N,FARRAY<int> IPVT,FARRAY<double> DET,FARRAY<double> WORK,int JOB)
        {
            #ifdef USING_BLAS
            dgedi_((double*)&A(1),&LDA,&N,(int*)&IPVT(1),(double*)&DET(1),(double*)&WORK(1),&JOB);
            #else






































































            A.array(LDA,N);
            IPVT.array(N);
            WORK.array(N);
            DET.array(2);

            double T,TEN;
            int I,J,K,KB,KP1,L,NM1;




            if (JOB/10 != 0) {
                DET(1) = 1.0;
                DET(2) = 0.0;
                TEN = 10.0;
                for (I=1; I<=N; ++I) {
                    if (IPVT(I) != I) DET(1) = -DET(1);
                    DET(1) = A(I,I)*DET(1);

                    if (DET(1) == 0.0) break;
                    while (fabs(DET(1)) < 1.00) {
                        DET(1) = TEN*DET(1);
                        DET(2) = DET(2) - 1.00;
                    }

                    while (fabs(DET(1)) >= TEN) {
                        DET(1) = DET(1)/TEN;
                        DET(2) = DET(2) + 1.00;
                    }
                }
            }



            if (JOB%10 == 0) return;
            for (K=1; K<=N; ++K) {
                A(K,K) = 1.00/A(K,K);
                T = -A(K,K);
                DSCAL(K-1,T,A(1,K),1);
                KP1 = K + 1;
                if (N >= KP1) {
                    for (J=KP1; J<=N; ++J) {
                        T = A(K,J);
                        A(K,J) = 0.0;
                        DAXPY(K,T,A(1,K),1,A(1,J),1);
                    }
                }
            }



            NM1 = N - 1;
            if (NM1 < 1) return;
            for (KB=1; KB<=NM1; ++KB) {
                K = N - KB;
                KP1 = K + 1;
                for (I=KP1; I<=N; ++I) {
                    WORK(I) = A(I,K);
                    A(I,K) = 0.0;
                }
                for (J=KP1; J<=N; ++J) {
                    T = WORK(J);
                    DAXPY(N,T,A(1,J),1,A(1,K),1);
                }
                L = IPVT(K);
                if (L != K) DSWAP(N,A(1,K),1,A(1,L),1);
            }
            return;
            #endif
        }

        void DSWAP(int N,FARRAY<double> DX,int INCX,FARRAY<double> DY,int INCY)
        {
            #ifdef USING_BLAS
            dswap_(&N,(double*)&DX(1),&INCX,(double*)&DY(1),&INCY);
            #else








































            DX.array(N);
            DY.array(N);
            double DTEMP1,DTEMP2,DTEMP3;

            if (N<=0) return;
            if (INCX!=INCY || INCX-1<0) {



                int IX = 1;
                int IY = 1;
                if (INCX<0) IX = (-N+1)*INCX + 1;
                if (INCY<0) IY = (-N+1)*INCY + 1;
                for (int I=1; I<=N; ++I) {
                    DTEMP1 = DX(IX);
                    DX(IX) = DY(IY);
                    DY(IY) = DTEMP1;
                    IX = IX + INCX;
                    IY = IY + INCY;
                }
                return;
            } else if (INCX-1==0) {






                int M = N%3;
                if ( M != 0 ) {
                    for (int I=1; I<=M; ++I) {
                        DTEMP1 = DX(I);
                        DX(I) = DY(I);
                        DY(I) = DTEMP1;
                    }
                    if ( N < 3 ) return;
                }
                int MP1 = M + 1;
                for (int I=MP1; I<=N; I+=3) {
                    DTEMP1 = DX(I);
                    DTEMP2 = DX(I+1);
                    DTEMP3 = DX(I+2);
                    DX(I) = DY(I);
                    DX(I+1) = DY(I+1);
                    DX(I+2) = DY(I+2);
                    DY(I) = DTEMP1;
                    DY(I+1) = DTEMP2;
                    DY(I+2) = DTEMP3;
                }
                return;
            } else {



                int NS = N*INCX;
                for (int I=1; I<=NS; I+=INCX) {
                    DTEMP1 = DX(I);
                    DX(I) = DY(I);
                    DY(I) = DTEMP1;
                }
            }
            #endif
        }

        double DDOT(int N,DFARRAY DX,int INCX,DFARRAY DY,int INCY)
        {
            #ifdef USING_BLAS
            return ddot_(&N,&DX(1),&INCX,&DY(1),&INCY);
            #else







































            DX.array(N);
            DY.array(N);

            double _DDOT = 0.0;
            if (N<=0) return _DDOT;
            if (INCX!=INCY || INCX<1) {



                int IX = 1;
                int IY = 1;
                if (INCX<0) IX = (-N+1)*INCX + 1;
                if (INCY<0) IY = (-N+1)*INCY + 1;
                for (int I=1; I<=N; ++I) {
                    _DDOT = _DDOT + DX(IX)*DY(IY);
                    IX = IX + INCX;
                    IY = IY + INCY;
                }
                return _DDOT;
            }
            if (INCX==1) {







                int M = N%5;
                if ( M != 0 ) {
                    for (int I=1; I<=M; ++I) {
                        _DDOT = _DDOT + DX(I)*DY(I);
                    }
                    if (N < 5 ) return _DDOT;
                }
                int MP1 = M + 1;
                for (int I=MP1; I<=N; I+=5) {
                    _DDOT = _DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
                            DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4);
                }
                return _DDOT;
            } else {



                int NS = N*INCX;
                for (int I=1; I<=NS; I+=INCX) {
                    _DDOT = _DDOT + DX(I)*DY(I);
                }

                return _DDOT;
            }
            #endif
        }

        void DGESL(DFARRAY A,int LDA,int N,IFARRAY IPVT,DFARRAY B,int JOB)
        {
            #ifdef USING_BLAS
            dgesl_(&A(1),&LDA,&N,&IPVT(1),&B(1),&JOB);
            #else




































































            IPVT.array(N);
            A.array(LDA,N);
            B.array(N);
            double T;
            int K,KB,L,NM1;

            NM1 = N - 1;
            if (JOB == 0) {




                if (NM1 >= 1) {
                    for (K=1; K<=NM1; ++K) {
                        L = IPVT(K);
                        T = B(L);
                        if (L != K) {
                            B(L) = B(K);
                            B(K) = T;
                        }
                        DAXPY(N-K,T,A(K+1,K),1,B(K+1),1);
                    }
                }



                for (KB=1; KB<=N; ++KB) {
                    K = N + 1 - KB;
                    B(K) = B(K)/A(K,K);
                    T = -B(K);
                    DAXPY(K-1,T,A(1,K),1,B(1),1);
                }
                return;
            }




            for (K=1; K<=N; ++K) {
                T = DDOT(K-1,A(1,K),1,B(1),1);
                B(K) = (B(K) - T)/A(K,K);
            }



            if (NM1 < 1) return;
            for (KB=1; KB<=NM1; ++KB) {
                K = N - KB;
                B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1);
                L = IPVT(K);
                if (L != K) {
                    T = B(L);
                    B(L) = B(K);
                    B(K) = T;
                }
            }
            return;
            #endif
        }


    }


}
