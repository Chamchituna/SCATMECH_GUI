











#include "random.h"
#include <limits>
#include "scatmech.h"
#include "matrixmath.h"
#include <iostream>
#include "fft.h"

#if _MSC_VER > 1
#include <process.h>
#else
#include <sys/types.h>
#include <unistd.h>
#endif

using namespace std;

namespace SCATMECH {



































































	static mutex_t mutex;

    namespace CMLIB {

        UNIrand::UNIrand()
        {
            M[0]=30788;
            M[1]=23052;
            M[2]=2053;
            M[3]=19346;
            M[4]=10646;
            M[5]=19427;
            M[6]=23975;
            M[7]=19049;
            M[8]=10949;
            M[9]=19693;
            M[10]=29746;
            M[11]=26748;
            M[12]=2796;
            M[13]=23890;
            M[14]=29168;
            M[15]=31924;
            M[16]=16499;
            I=5;
            J=17;
            M1=32767;
            M2=256;
            seeded=false;
        }

        void UNIrand::seed(int JD)
        {
            int MDIG = numeric_limits<int>::digits+1;
            if (MDIG < 16) throw SCATMECH_exception("UNI--MDIG LESS THAN 16");
            M1= (1<<(MDIG-2)) + ((1<<(MDIG-2))-1);
            M2 = (1<<(MDIG/2));
            int iabsJD = abs(JD);
            int JSEED = iabsJD<M1 ? iabsJD : M1;
            if ( JSEED%2 == 0 ) JSEED=JSEED-1;
            int K0 = 9069%M2;
            int K1 = 9069/M2;
            int J0 = JSEED%M2;
            int J1 = JSEED/M2;
            for (I=1; I<=17; ++I) {
                JSEED = J0*K0;
                J1 = (JSEED/M2+J0*K1+J1*K0) % (M2/2);
                J0 = JSEED % M2;
                M[I-1] = J0+M2*J1;
            }
            I=5;
            J=17;
            seeded=true;
        }

        double UNIrand::operator()()
        {
			if (!seeded) {
				mutex.lock();
				seed((int)time(NULL) + (seed_add += seed_adder));
				mutex.unlock();
			}
            int K=M[I-1]-M[J-1];
            if (K < 0) K=K+M1;
            M[J-1]=K;
            I=I-1;
            if (I == 0) I=17;
            J=J-1;
            if (J == 0) J=17;
            double result = (double)K/(double)M1;
            return result;
        }
    }

    Random_Number::Random_Number()
    {
        initialize();
        *this = uniform(0.,1.);
    }

    Random_Number::Random_Number(const SCATMECH::Table& distrib) {
        initialize();
        set_distribution(distrib);
    }

    Table Random_Number::uniform(double xmin,double xmax) {
        Table::VECTOR distrib;
        distrib.push_back(std::pair<double,double>(xmin,1.));
        distrib.push_back(std::pair<double,double>(xmax,1.));
        return Table(distrib);
    }

    Table Random_Number::normal(double stdev, int npoints, double nstdev) {
        Table::VECTOR normal;
        for (int i=0; i<npoints; ++i) {
            double x = (i-npoints/2.)/(double)(npoints)*2.*nstdev*stdev;
            double xnorm = x/stdev;
            double y = exp(-xnorm*xnorm/2.);
            normal.push_back(std::pair<double,double>(x,y));
        }
        return Table(normal);
    }





















    Table Random_Number::sinesqr(int npoints) {
        Table::VECTOR prob;
        for (int i=0; i<npoints; ++i) {
            double x = i/(double)(npoints)*pi;
            double y = sqr(sin(x));
            prob.push_back(std::pair<double,double>(x,y));
        }
        return Table(prob);
    }

    Table Random_Number::exponential(double decay, int npoints,double ndecay) {
        SCATMECH::Table::VECTOR exponential;
        for (int i=0; i<npoints; ++i) {
            double x = (double)i/(double)(npoints)*decay*ndecay;
            double xnorm = x/decay;
            double y = exp(-xnorm);
            exponential.push_back(std::pair<double,double>(x,y));
        }
        return Table(exponential);
    }

    Table Random_Number::log_normal(double mode, double stdev, int npoints, double nstdev) {
        SCATMECH::Table::VECTOR log_normal;
        double log10mode = log10(mode);
        for (int i=0; i<npoints; ++i) {
            double log10x =(i-npoints/2.)/(double)(npoints)*2.*nstdev*stdev+log10mode;
            double x = pow(10.,log10x);
            double xnorm = (log10x-log10mode)/stdev;
            double y = exp(-xnorm*xnorm/2.)/x;
            log_normal.push_back(std::pair<double,double>(x,y));
        }
        return Table(log_normal);
    }

    double Random_Number::operator()() {
        return lookup.value(rand());
    }

    double Random_Number::operator()(double offset) {
        double r = rand();
        if (r<offset) {
            Table::VECTOR table = lookup.get_table();
            double a = table[0].second;
            double b = table[table.size()-1].second;
            return r/offset*(b-a)+a;
        } else {
            return lookup.value((r-offset)/(1.-offset));
        }
    }

    void Random_Number::set_seed(long seed) {
        rand.seed(seed);
    }

    void Random_Number::set_distribution(const SCATMECH::Table& distrib) {
        Table::VECTOR temp;

        double integral=0;
        double x=distrib.get_table().begin()->first;
        double y=distrib.get_table().begin()->second;

        temp.push_back(std::pair<double,double>(0.,x));

        SCATMECH::Table::VECTOR::const_iterator cp;
        for (cp=distrib.get_table().begin(),++cp; cp<distrib.get_table().end(); ++cp) {
            double dx = cp->first - x;
            integral += (cp->second + y)/2.*dx;

            x = cp->first;
            y = cp->second;

            temp.push_back(std::pair<double,double>(integral,x));
        }

        SCATMECH::Table::VECTOR::iterator p;
        for (p=temp.begin(); p<temp.end(); ++p) {
            p->first /= integral;
        }
        lookup = Table(temp);
    }

    void Random_Number::initialize() {
		mutex.lock();
        set_seed(-(int)time(NULL)+(seed_add+=seed_adder));
		mutex.unlock();
    }

    int Random_Number::seed_add=12345;
    int CMLIB::UNIrand::seed_add=123456;

    namespace {
        int mygetpid() {
            #if _MSC_VER > 1000
            int result =_getpid();
            #else
            int result = getpid();
            #endif
            return result;
        }
    }

    int Random_Number::seed_adder=594*mygetpid();
    int CMLIB::UNIrand::seed_adder=594*mygetpid();

    Multivariate_Random_Number::Multivariate_Random_Number()
    {
        random = Random_Number::normal();
        N=0;
    }

    Multivariate_Random_Number::Multivariate_Random_Number(int n, const std::vector<double>& covmatrix)
    {
        N=n;
        random = Random_Number::normal();
        if (covmatrix.size()!=n*n) throw SCATMECH_exception("Invalid size for covariance matrix");

        vector<COMPLEX> cov(n*n);

        for (int i=0; i<n; ++i) {
            for (int j=0; j<n; ++j) {
                cov[i+n*j] = covmatrix[i+n*j];
            }
        }
        eigenvalues.resize(N);
        eigenvectors.resize(N*N);
        eigen(cov,eigenvalues,eigenvectors,N);
        temp.resize(N);
        for (int k=0; k<n; ++k) eigenvalues[k] = sqrt(eigenvalues[k]);
    }

    void Multivariate_Random_Number::get(std::vector<double>& values)
    {
        values.resize(N);
        for (int i=0; i<N; ++i) {
            temp[i] = random()*real(eigenvalues[i]);
        }
        for (int j=0; j<N; ++j) {
            values[j]=0;
            for (int k=0; k<N; ++k) {
                values[j] += real(eigenvectors[j+N*k])*temp[k];
            }
        }
    }

    Table Random_Curve(double rms,double corrlength, double exponent, double period, int n, int seed)
    {
        static Random_Number rand = Random_Number::uniform(0.,2*pi);
        if (seed!=0) rand.set_seed(seed);

        int N = n;
        valarray<COMPLEX> fC(N);
        vector<double> f(N);
        vector<double> fx(N);

        int i;
        double dx = period/N;




        for (i=0; i<N; ++i) {
            double z = (i<N/2) ? i*dx : (N-i)*dx;
            fC[i]=exp(-pow(sqr(z/corrlength),exponent));
        }


        fft1d(fC,+1);


        for (i=0; i<N; ++i) fC[i]=sqrt(fC[i]);



        COMPLEX phase = exp(COMPLEX(0,1)*rand());
        fC[N/2]*=phase;
        fC[0]*=conj(phase);

        for (i=1; i<N/2; ++i) {
            phase = exp(COMPLEX(0,1)*rand());
            fC[i] *= phase;
            fC[N-i] *= conj(phase);
        }


        fft1d(fC,-1);


        double sum=0;
        double sum2=0;
        for (i=0; i<N; ++i) {
            f[i]=real(fC[i])/N;
            fx[i]=i*dx;
            sum+=f[i];
            sum2+=sqr(f[i]);
        }
        double mean=sum/N;
        double r=sqrt(sum2/N-sqr(mean));
        for (i=0; i<N; ++i) {
            f[i]=(f[i]-mean)*rms/r;
        }


        return Table(fx,f);
    }

}
