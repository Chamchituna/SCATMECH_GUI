











#ifndef SCATMECH_REFLECTANCE_H
#define SCATMECH_REFLECTANCE_H

#include "inherit.h"
#include "scattabl.h"
#include <string>

namespace SCATMECH {





    class Reflectance : public Model {
        public:




            virtual double Get_Reflectance(double _lambda) {
                lambda = _lambda;
                if (lambda!=old_lambda) {
                    set_recalc(1);
                    old_lambda=lambda;
                }
                SETUP();
                return reflectance;
            }
        protected:


            double lambda;


            double old_lambda;



            double reflectance;

            DECLARE_MODEL();
    };

    typedef Model_Ptr<Reflectance> Reflectance_Ptr;
    void Register(const Reflectance* x);








    class Kubelka_Munk_Reflectance : public Reflectance {
        protected:

            void setup();

            DECLARE_MODEL();


            DECLARE_PARAMETER(double,thickness);


            DECLARE_PARAMETER(Table,absorption);


            DECLARE_PARAMETER(Table,scattering);


            DECLARE_PARAMETER(Reflectance_Ptr,substrate);
    };





    class Table_Reflectance : public Reflectance {
        protected:

            void setup();

            DECLARE_MODEL();


            DECLARE_PARAMETER(Table,table);
    };




    class Equation_Reflectance : public Reflectance {
        protected:

            void setup();

            DECLARE_MODEL();



            DECLARE_PARAMETER(std::string,expression);
    };
}

#endif
