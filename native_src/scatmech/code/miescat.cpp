












#include "scatmech.h"
#include "miescat.h"

















using namespace std;


namespace SCATMECH {


    MieScatterer::
    MieScatterer()
    {
        NSTOP=0;
        NMX=0;
    }

    void
    MieScatterer::
    setup()
    {
        SphericalScatterer::setup();


        COMPLEX REFREL = m;

        COMPLEX Y=x*REFREL;
        double ymod=abs(Y);



        double XSTOP=x+4.*pow(x,0.3333)+2.;

        NMX=(int)(((XSTOP>ymod)? XSTOP : ymod )+15);


        NSTOP=(int)(XSTOP);





        A.resize(NMX);
        B.resize(NMX);
        E.resize(NMX);
        F.resize(NMX);



        vector<COMPLEX> D;
        D.resize(NMX);
        D[NMX-1]=COMPLEX(0.,0.);
        int NN=NMX-1;
        int N;

        for (N=0; N<NN; ++N) {
            E[N]=NMX-N;
            D[NMX-N-2]=(E[N]/Y)-(1./(D[NMX-N-1]+E[N]/Y));
        }



        double  PSI0=cos(x);
        double  PSI1=sin(x);
        double  CHI0=-sin(x);
        double  CHI1=cos(x);
        COMPLEX XI1=COMPLEX(PSI1,-CHI1);

        for (N=0; N<NSTOP; ++N) {
            E[N]=N+1;
            F[N]=(2.*E[N]+1.)/(E[N]*(E[N]+1.));




            double PSI=(2.*E[N]-1.)*PSI1/x-PSI0;
            double CHI=(2.*E[N]-1.)*CHI1/x-CHI0;
            COMPLEX XI=COMPLEX(PSI,-CHI);

            A[N]=(D[N]/REFREL+E[N]/x)*PSI-PSI1;
            A[N]=A[N]/((D[N]/REFREL+E[N]/x)*XI-XI1);
            B[N]=(REFREL*D[N]+E[N]/x)*PSI-PSI1;
            B[N]=B[N]/((REFREL*D[N]+E[N]/x)*XI-XI1);

            PSI0=PSI1;
            PSI1=PSI;
            CHI0=CHI1;
            CHI1=CHI;
            XI1=COMPLEX(PSI1,-CHI1);
        }
    }




    COMPLEX
    MieScatterer::
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
    MieScatterer::
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
    MieScatterer::
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
    MieScatterer::
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
    MieScatterer::
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

    DEFINE_MODEL(MieScatterer,SphericalScatterer,
                 "Scattering by a homogeneous sphere.");



}



