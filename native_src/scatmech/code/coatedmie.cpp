












#include "scatmech.h"
#include "coatedmie.h"
#include "bobvlieg.h"
#include "askuser.h"

using namespace std;


namespace SCATMECH {

    using namespace BobVlieg_Supp;




    CoatedMieScatterer::
    CoatedMieScatterer()
    {
        NSTOP=0;
        NMX=0;
    }

    void
    CoatedMieScatterer::
    setup()
    {
        SphericalScatterer::setup();

        COMPLEX REFREL = m;

        COMPLEX Y=x*REFREL;
        double ymod=abs(Y);



        double XSTOP=x+4.*pow(x,0.3333)+2.;

        if (nmax<0) {
            XSTOP = XSTOP + nmax;
        }
        if (nmax>0) {
            XSTOP = nmax;
        }

        NMX=(int)(((XSTOP>ymod)? XSTOP : ymod )+15);

        NSTOP=(int)(XSTOP);

        A.resize(NMX);
        B.resize(NMX);
        E.resize(NMX);
        F.resize(NMX);


        int N;

        double x = k*(radius-thickness);
        double y = k*(radius);

        COMPLEX m2 = (COMPLEX)coating.index(lambda)/(COMPLEX)medium.index(lambda);
        COMPLEX m1 = (COMPLEX)sphere.index(lambda)/(COMPLEX)medium.index(lambda);




        for (N=NSTOP-1; N>=0; --N) {
            E[N]=N+1;
            F[N]=(2.*E[N]+1.)/(E[N]*(E[N]+1.));

            int n=N+1;
            COMPLEX An = (m2*psi(n,m2*x)*psi_(n,m1*x)-
                          m1*psi_(n,m2*x)*psi(n,m1*x))/
                         (m2*chi(n,m2*x)*psi_(n,m1*x)-
                          m1*chi_(n,m2*x)*psi(n,m1*x));
            COMPLEX an = (psi(n,y)*(psi_(n,m2*y)-An*chi_(n,m2*y)) -
                          m2*psi_(n,y)*(psi(n,m2*y)-An*chi(n,m2*y)))/
                         (zeta(n,y)*(psi_(n,m2*y)-An*chi_(n,m2*y)) -
                          m2*zeta_(n,y)*(psi(n,m2*y)-An*chi(n,m2*y)));

            COMPLEX Bn = (m2*psi(n,m1*x)*psi_(n,m2*x)-
                          m1*psi(n,m2*x)*psi_(n,m1*x))/
                         (m2*chi_(n,m2*x)*psi(n,m1*x)-
                          m1*psi_(n,m1*x)*chi(n,m2*x));
            COMPLEX bn = (m2*psi(n,y)*(psi_(n,m2*y)-Bn*chi_(n,m2*y))-
                          psi_(n,y)*(psi(n,m2*y)-Bn*chi(n,m2*y)))/
                         (m2*zeta(n,y)*(psi_(n,m2*y)-Bn*chi_(n,m2*y))-
                          zeta_(n,y)*(psi(n,m2*y)-Bn*chi(n,m2*y)));

            A[N]=an;
            B[N]=bn;
        }
    }




    COMPLEX
    CoatedMieScatterer::
    s1(double theta)
    {
        SETUP();

        COMPLEX S1=0.;
        double PI0=0.;
        double PI1=1.;
        double AMU=cos(theta);
        int N;
        for (N=0; N<NSTOP; ++N) {
            double PPI=PI1;
            double TAU=E[N]*AMU*PPI-(E[N]+1.)*PI0;
            S1+=F[N]*(A[N]*PPI+B[N]*TAU);
            PI1=((2.*E[N]+1.)*AMU*PPI-(E[N]+1.)*PI0)/E[N];
            PI0=PPI;
        }
        return S1;
    }




    COMPLEX
    CoatedMieScatterer::
    s2(double theta)
    {
        SETUP();

        COMPLEX S2=0.;
        double PI0=0.;
        double PI1=1.;
        double AMU=cos(theta);
        int N;
        for (N=0; N<NSTOP; ++N) {
            double PPI=PI1;
            double TAU=E[N]*AMU*PPI-(E[N]+1.)*PI0;
            S2+=F[N]*(A[N]*TAU+B[N]*PPI);
            PI1=((2.*E[N]+1.)*AMU*PPI-(E[N]+1.)*PI0)/E[N];
            PI0=PPI;
        }
        return S2;
    }




    double
    CoatedMieScatterer::
    Csca()
    {
        SETUP();

        double result=0;
        for (int i=0; i<NSTOP; ++i) {
            int n=i+1;

            result += (2.*n+1.)*(norm(A[i])+norm(B[i]));
        }
        result *= 2.*pi/sqr(k);
        return result;
    }




    double
    CoatedMieScatterer::
    Cext()
    {
        SETUP();

        double result=0;
        for (int i=0; i<NSTOP; ++i) {
            int n=i+1;

            result += (2.*n+1.)*real(A[i]+B[i]);
        }
        result *= 2.*pi/sqr(k);
        return result;
    }




    double
    CoatedMieScatterer::
    Cback()
    {
        SETUP();

        COMPLEX result1=0;
        for (int i=0; i<NSTOP; ++i) {
            int n=i+1;

            result1 += ((2.*n+1.)*pow(-1.,(double)n))*(A[i]-B[i]);
        }
        double result = pi/sqr(k)*norm(result1);
        return result;
    }

    DEFINE_MODEL(CoatedMieScatterer,SphericalScatterer,
                 "Scattering from a homogeneous sphere with a coating.");

    DEFINE_PARAMETER(CoatedMieScatterer,dielectric_function,coating,"Sphere coating material","(1,0)",0xFF);

    DEFINE_PARAMETER(CoatedMieScatterer,double,thickness,"Sphere coating thickness [um]","0",0xFF);
    DEFINE_PARAMETER(CoatedMieScatterer,int,nmax,"Highest spherical harmonic order (0 = use BH value)","0",0xFF);


}



