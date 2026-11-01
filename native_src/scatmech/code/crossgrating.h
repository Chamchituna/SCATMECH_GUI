











#ifndef CROSSGRATING_H
#define CROSSGRATING_H

#include <vector>
#include <algorithm>
#include "matrixmath.h"
#include "inherit.h"
#include "dielfunc.h"
#include "grating.h"
#include "scattabl.h"

namespace SCATMECH {

    class CrossGrating : public Model {
        public:

            CrossGrating() : order1(1),order2(1),lambda(0.5) {}






			CFARRAY get_EPS0() {
                SETUP();
                return EPS0;
            }

            CFARRAY get_EPS11() {
                SETUP();
                return EPS11;
            }

			CFARRAY get_EPS12() {
				SETUP();
				return EPS12;
			}


            CFARRAY get_EPS2() {
                SETUP();
                return EPS2;
            }

            CFARRAY get_EPS3() {
                SETUP();
                return EPS3;
            }


			CFARRAY get_MU0() {
				SETUP();
				return MU0;
			}

			CFARRAY get_MU11() {
				SETUP();
				return MU11;
			}

			CFARRAY get_MU12() {
				SETUP();
				return MU12;
			}


			CFARRAY get_MU2() {
				SETUP();
				return MU2;
			}

			CFARRAY get_MU3() {
				SETUP();
				return MU3;
			}



            double get_thick(int level) {
                SETUP();
                return thick(level);
            }


            int get_levels() {
                SETUP();
                return levels;
            }


            void set_lambda(double _lambda) {
                lambda = _lambda;
                set_recalc(0xFF);
            }
            double get_lambda() const {
                return lambda;
            }


            void set_order1(int _order1) {
                order1 = _order1;
                set_recalc(0xFF);
            }
            int get_order1() {
                return order1;
            }


            void set_order2(int _order2) {
                order2 = _order2;
                set_recalc(0xFF);
            }
            int get_order2() {
                return order2;
            }


            double get_zeta() {
                SETUP();
                return zeta;
            }


            double get_d1() {
                SETUP();
                return d1;
            }


            double get_d2() {
                SETUP();
                return d2;
            }

        protected:

            void setup();


            virtual void set_parameter_base(
                const STRING& parameter,
                const STRING& value
            );



            CFARRAY EPS0;
            CFARRAY EPS11;
			CFARRAY EPS12;
			CFARRAY EPS2;
            CFARRAY EPS3;
			CFARRAY MU0;
			CFARRAY MU11;
			CFARRAY MU12;
			CFARRAY MU2;
			CFARRAY MU3;


            int levels;
            double zeta;
            double d1;
            double d2;
            DFARRAY thick;


            double lambda;
            int order1;
            int order2;

        protected:

            DECLARE_MODEL();
            DECLARE_PARAMETER(dielectric_function,medium_i);
            DECLARE_PARAMETER(dielectric_function,medium_t);
    };





    class Gridded_CrossGrating : public CrossGrating {
        public:
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,zeta);
            DECLARE_PARAMETER(double,d1);
            DECLARE_PARAMETER(double,d2);
            DECLARE_PARAMETER(int,grid1);
            DECLARE_PARAMETER(int,grid2);

		public:
			Gridded_CrossGrating() : nonmagnetic(true), isotropic(true) {}
        protected:
            void setup() {
                CrossGrating::setup();

                CrossGrating::zeta = zeta;
                CrossGrating::d1 = d1;
                CrossGrating::d2 = d2;

                coszeta=cos(zeta*deg);
                sinzeta=sin(zeta*deg);
            }


            virtual void set_parameter_base(
                const STRING& parameter,
                const STRING& value
            );



            CFARRAY eps;



			CFARRAY eps1, eps2, eps3;
			CFARRAY mu;
			CFARRAY mu1, mu2, mu3;


            void FourierFactorize();


			void FourierFactorize(CFARRAY& eps1, CFARRAY& eps2, CFARRAY& eps3, CFARRAY& _EPS0, CFARRAY& _EPS11, CFARRAY& _EPS12, CFARRAY& _EPS2, CFARRAY& _EPS3);


            COMPLEX epsilon(const dielectric_function& e) {
                return (COMPLEX)(e.epsilon(lambda));
            }

            double coszeta,sinzeta;

            void getxy(int i,int j,double &x, double &y);

			bool isotropic;
			bool nonmagnetic;
    };

    void Register(const CrossGrating* x);
    typedef Model_Ptr<CrossGrating> CrossGrating_Ptr;

}

#endif
