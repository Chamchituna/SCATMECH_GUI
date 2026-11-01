











#ifndef TMATRIX_H
#define TMATRIX_H

#include "inherit.h"
#include "sphrscat.h"
#include "scatmatrix.h"

namespace SCATMECH {

    class Axisymmetric_Shape;
    typedef Model_Ptr<Axisymmetric_Shape> Axisymmetric_Shape_Ptr;
    void Register(const Axisymmetric_Shape* x);

    class shape_rec {
        public:
            double shape;
            double dshape;
            double abs;
            double wt;
    };

    typedef std::vector<shape_rec> shape_vec;

    class Axisymmetric_Shape : public Model {
        public:

            double Get_Volume();



            double Get_MaxRadius();



            double Get_Base_Length();



            double Get_Breast();


            void Get_TMatrix(ScatterTMatrix& T, int lmax, int mmax, double k, COMPLEX& index);


            void Write(const std::string& filename);

        protected:
            void Insert_Points(double start,double end);
            void Flip();
            void Renormalize(double newV);


            virtual void set_parameter_base(const STRING& parameter, const STRING& value);

            DECLARE_MODEL();
            DECLARE_PARAMETER(int,npoints);

        private:
            shape_vec shape;

        protected:
            virtual void setup();
            virtual double f(double theta) const = 0;
            virtual double df(double theta) const;
    };

    class Ellipsoid_Axisymmetric_Shape : public Axisymmetric_Shape {
            void setup();
            double f(double theta) const;
            double df(double theta) const;

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,vertical);
            DECLARE_PARAMETER(double,horizontal);
            DECLARE_PARAMETER(double,offset);
    };

    class Cylinder_Axisymmetric_Shape : public Axisymmetric_Shape {
            void setup();
            double f(double theta) const;
            double df(double theta) const;

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,radius);
            DECLARE_PARAMETER(double,length);
            DECLARE_PARAMETER(double,corner);
            DECLARE_PARAMETER(int,renorm);
        private:
            double theta0,theta1,theta2,c;
    };

    class Table_Axisymmetric_Shape : public Axisymmetric_Shape {
            void setup();
            double f(double theta) const;

            DECLARE_MODEL();
            DECLARE_PARAMETER(Table,table);
            DECLARE_PARAMETER(double,scale);
    };

    class Conical_Axisymmetric_Shape : public Axisymmetric_Shape {
            void setup();
            double f(double theta) const;
            double df(double theta) const;

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,radius);
            DECLARE_PARAMETER(double,height);
            DECLARE_PARAMETER(double,offset);
            DECLARE_PARAMETER(double,corner_base);
            DECLARE_PARAMETER(double,corner_apex);
            DECLARE_PARAMETER(int,updown);
            DECLARE_PARAMETER(int,renorm);
        private:
            double z0,beta,a0,theta0,theta1,theta2,theta3,c0,g;
    };

    class Indented_Ellipsoid_Axisymmetric_Shape : public Axisymmetric_Shape {
            void setup();
            double f(double theta) const;

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,vertical);
            DECLARE_PARAMETER(double,horizontal);
            DECLARE_PARAMETER(double,indent);
            DECLARE_PARAMETER(int,updown);
            DECLARE_PARAMETER(int,renorm);
        private:
            double theta0,B;
    };

    class Double_Conical_Axisymmetric_Shape : public Axisymmetric_Shape {
            void setup();
            double f(double theta) const;

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,height);
            DECLARE_PARAMETER(double,radius);
            DECLARE_PARAMETER(double,corner_apex);
            DECLARE_PARAMETER(double,corner_waist);
            DECLARE_PARAMETER(int,renorm);
        private:
            double theta0,theta1,alpha,beta,c1,c2;
    };

    class Chebyshev_Axisymmetric_Shape : public Axisymmetric_Shape {
            void setup();
            double f(double theta) const;

            DECLARE_MODEL();
            DECLARE_PARAMETER(int,n);
            DECLARE_PARAMETER(double,start);
            DECLARE_PARAMETER(double,end);
            DECLARE_PARAMETER(double,amplitude);
            DECLARE_PARAMETER(double,vertical);
            DECLARE_PARAMETER(double,horizontal);
            DECLARE_PARAMETER(int,renorm);
        private:
            double a;
    };

    class Gauss_Legendre_Integration
    {
        public:
            Gauss_Legendre_Integration(int order,double x1,double x2);
            double weight(int i) const {
                return w[i];
            }
            double abscissa(int i) const {
                return x[i];
            }
            int n() const {
                return x.size();
            }
        private:
            std::vector<double> w;
            std::vector<double> x;
    };


}
#endif
