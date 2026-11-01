











#ifndef RANDOM_H
#define RANDOM_H

#include <ctime>
#include <vector>
#include "scattabl.h"
#include <cmath>

namespace SCATMECH {

    namespace CMLIB {




        class UNIrand {
            public:
                UNIrand();
                void seed(int JD);
                double operator()();

            private:
                static int seed_add;
                static int seed_adder;
                int M[17],I,J,M1,M2;
                bool seeded;
        };
    }



    class Random_Number {
        public:



            Random_Number();




            Random_Number(const SCATMECH::Table& distrib);






            static Table uniform(double xmin=0.,double xmax=1.);


            static Table normal(double stdev = 1., int npoints=400, double nstdev=5.);


            static Table sinesqr(int npoints=400);


            static Table exponential(double decay, int npoints=400,double ndecay=10.);



            static Table log_normal(double mode, double stdev = 1., int npoints=400, double nstdev=5.);




            double operator()();





            double operator()(double offset);




            void set_seed(long seed);

        private:

            void set_distribution(const SCATMECH::Table& pdf);
            void initialize();
            SCATMECH::Table lookup;
            CMLIB::UNIrand rand;
            static int seed_add;
            static int seed_adder;
    };





    class Multivariate_Random_Number {
        public:



            Multivariate_Random_Number();




            Multivariate_Random_Number(int n, const std::vector<double>& covmatrix);




            void get(std::vector<double>& values);




            void set_seed(long seed) {
                random.set_seed(seed);
            }

        private:
            Random_Number random;
            std::vector<COMPLEX> eigenvalues,eigenvectors;
            std::vector<double> temp;
            int N;
    };





    Table Random_Curve(double rms,
                       double corrlength,
                       double exponent,
                       double period,
                       int n,
                       int seed=0
                      );
}

#endif
