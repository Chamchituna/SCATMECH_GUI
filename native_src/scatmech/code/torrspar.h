











#ifndef SCATMECH_TORRSPAR_H
#define SCATMECH_TORRSPAR_H

#include "scatmech.h"
#include "facet.h"


namespace SCATMECH {


    class Shadow_Function : public Model {
        public:
            virtual double f(double thetai,double thetas,double phis) = 0;
            DECLARE_MODEL();
    };

    class Unit_Shadow_Function: public Shadow_Function
    {
        public:
            double f(double thetai,double thetas,double phis) {
                return 1.;
            }
            DECLARE_MODEL();
    };

    class Torrance_Sparrow_Shadow_Function: public Shadow_Function
    {
        public:
            double f(double thetai,double thetas,double phis);
            DECLARE_MODEL();
    };






    class Smith_Shadow_Function: public Shadow_Function
    {
        public:
            double f(double thetai,double thetas,double phis);

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,w);

        private:
            double S(double theta);
    };









    class Maxwell_Beard_Shadow_Function: public Shadow_Function
    {
        public:
            double f(double thetai,double thetas,double phis);

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,tau);
            DECLARE_PARAMETER(double,Omega);
    };

    class Table_Shadow_Function: public Shadow_Function
    {
        public:
            double f(double thetai,double thetas,double phis);
            DECLARE_MODEL();
            DECLARE_PARAMETER(Table,T);
    };

    typedef Model_Ptr<Shadow_Function> Shadow_Function_Ptr;







    class Shadowed_Facet_BRDF_Model: public Facet_BRDF_Model
    {
        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(Shadow_Function_Ptr,shadow);

        protected:

            JonesMatrix jones();
    };

    void Register(const Shadow_Function* x);


}


#endif
