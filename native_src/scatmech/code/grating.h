











#ifndef SCATMECH_GRATING_H
#define SCATMECH_GRATING_H

#include "inherit.h"
#include "dielfunc.h"
#include "filmtran.h"
#include <map>

namespace SCATMECH {

    class RCW_Model;

    class Grating : public Model {
        public:
            Grating() {
                anisotropic = false;
                magnetic = false;
                lambda=1.;
            }

            virtual COMPLEX fourier(int order,int level,int recip=0)
            {
                SETUP();
                if (anisotropic) error("Calling fourier with an anisotropic grating");
                return fourierx(order,level,recip);
            }





            virtual COMPLEX fourierx(int order,int level,int recip=0);
            virtual COMPLEX fouriery(int order,int level,int recip=0);
            virtual COMPLEX fourierz(int order,int level,int recip=0);





            virtual COMPLEX fouriermux(int order,int level,int recip=0);
            virtual COMPLEX fouriermuy(int order,int level,int recip=0);
            virtual COMPLEX fouriermuz(int order,int level,int recip=0);


            int get_levels() {
                SETUP();
                return thickness.size();
            }


            double get_thickness(int level) {
                SETUP();
                return thickness[level];
            }



            virtual COMPLEX eps(double x,int level,int direction);



            COMPLEX mu(double x,int level,int direction);


            double get_total_thickness() {
                SETUP();
                double result=0;
                for (int i=0; i<(int)thickness.size(); ++i) result += thickness[i];
                return result;
            }






            int get_level(double z) {
                SETUP();
                if (z>0) return -1;
                for (int i=0; i<(int)thickness.size(); ++i) {
                    z += thickness[i];
                    if (z>0) return i;
                }
                return -2;
            }

            void set_lambda(double _lambda) {
                lambda = _lambda;
                set_recalc(0x01);
            }
            double get_lambda() const {
                return lambda;
            }

            bool is_anisotropic() {
                SETUP();
                return anisotropic;
            }
            bool is_magnetic() {
                SETUP();
                return magnetic;
            }

            friend class RCW_Model;

            const std::vector<std::vector<double> > & get_position() {
                SETUP();
                return position;
            }
            const std::vector<std::vector<COMPLEX> > & get_materialx() {
                SETUP();
                return materialx;
            }
            const std::vector<std::vector<COMPLEX> > & get_materialy() {
                SETUP();
                return materialy;
            }
            const std::vector<std::vector<COMPLEX> > & get_materialz() {
                SETUP();
                return materialz;
            }
            const std::vector<std::vector<COMPLEX> > & get_materialmux() {
                SETUP();
                return materialmux;
            }
            const std::vector<std::vector<COMPLEX> > & get_materialmuy() {
                SETUP();
                return materialmuy;
            }
            const std::vector<std::vector<COMPLEX> > & get_materialmuz() {
                SETUP();
                return materialmuz;
            }

            const std::vector<double> & get_thickness() {
                SETUP();
                return thickness;
            }

        protected:

            bool anisotropic;
            bool magnetic;

            typedef std::vector<std::vector<COMPLEX> > epsvector_t;


            std::vector<std::vector<double> > position;
            epsvector_t materialx,materialy,materialz;
            epsvector_t materialmux,materialmuy,materialmuz;
            std::vector<double> thickness;


            COMPLEX epsilon(const dielectric_function& e);

            double lambda;

        protected:

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,period);
            DECLARE_PARAMETER(dielectric_function,medium_i);
            DECLARE_PARAMETER(dielectric_function,medium_t);
    };

    void Register(const Grating* x);
    typedef Model_Ptr<Grating> Grating_Ptr;

    class Fourier_Component_Calculator {
        public:
            Fourier_Component_Calculator(double _period,int _order,int _recip);
            void enter(double x,COMPLEX y);
            COMPLEX result() {
                return sum;
            }

        private:
            double period,last_x;
            int order,recip;
            COMPLEX sum,last_y;
            bool begin;
    };

    class Generic_Grating : public Grating {
        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(std::string,filename);
            DECLARE_PARAMETER(std::string,pstring);
            DECLARE_PARAMETER(int,nlayers);

        public:

            struct segment {
                COMPLEX vertex1,vertex2;
                std::string mat1,mat2;
            };
            typedef std::vector<segment> segmentvector;

            struct bounds {
                double x;
                std::string mat1,mat2;

                bool operator>(const bounds& b) const {
                    return x>b.x;
                }
                bool operator<(const bounds& b) const {
                    return x<b.x;
                }
            };

            const segmentvector& Get_Boundaries() {
                SETUP();
                return boundaries;
            }

            void print_parameters(std::ostream& os,const std::string& prefix="") const;
            void print_ggparameters(std::ostream& os) const;
            void print_variables(std::ostream& os) const;
            void print_materials(std::ostream& os) const;
            void print_boundaries(std::ostream& os) const;

        protected:
            void setup();
            void error(const std::string& message) const;


            virtual STRING get_parameter_base(const STRING& parameter) const;
            virtual void set_parameter_base(const STRING& parameter, const STRING& value);

        private:

            typedef std::map<std::string,COMPLEX> epsmap_t;
            epsmap_t epsx,epsy,epsz;
            epsmap_t mux,muy,muz;

            std::string filecontents;
            std::string old_filename;

            typedef std::map<std::string,double> varsmap;
            varsmap vars;

            varsmap parameters;
            varsmap override;

            segmentvector boundaries;

            void gg_error(const std::string& message) const;

            COMPLEX get_complex_value(std::istream& is) const;
    };

    class Dielectric_Stack_Grating : public Grating {
        protected:
            void setup();

        private:
            DECLARE_MODEL();
            DECLARE_PARAMETER(StackModel_Ptr,stackepsx);
            DECLARE_PARAMETER(StackModel_Ptr,stackepsy);
            DECLARE_PARAMETER(StackModel_Ptr,stackepsz);
            DECLARE_PARAMETER(StackModel_Ptr,stackmux);
            DECLARE_PARAMETER(StackModel_Ptr,stackmuy);
            DECLARE_PARAMETER(StackModel_Ptr,stackmuz);
    };

    class Single_Line_Grating : public Grating {
        protected:
            void setup();

        private:
            DECLARE_MODEL();
            DECLARE_PARAMETER(dielectric_function,material);
            DECLARE_PARAMETER(dielectric_function,space);
            DECLARE_PARAMETER(double,topwidth);
            DECLARE_PARAMETER(double,bottomwidth);
            DECLARE_PARAMETER(double,offset);
            DECLARE_PARAMETER(int,nlevels);
            DECLARE_PARAMETER(double,height);
    };

    class Corner_Rounded_Grating : public Grating {
        protected:
            void setup();

        private:
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,width);
            DECLARE_PARAMETER(double,height);
            DECLARE_PARAMETER(dielectric_function,material);
            DECLARE_PARAMETER(double,sidewall);
            DECLARE_PARAMETER(double,radiusb);
            DECLARE_PARAMETER(double,radiust);
            DECLARE_PARAMETER(int,nlevels);

    };

    class Triangular_Grating : public Grating {
        protected:
            void setup();

        private:
            DECLARE_MODEL();
            DECLARE_PARAMETER(dielectric_function,material);
            DECLARE_PARAMETER(double,amplitude);
            DECLARE_PARAMETER(double,aspect);
            DECLARE_PARAMETER(int,nlevels);
    };

    class Sinusoidal_Relief_Grating : public Grating {
        protected:
            void setup();

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,amplitude);
            DECLARE_PARAMETER(double,base);
            DECLARE_PARAMETER(dielectric_function,material);
            DECLARE_PARAMETER(int,nlevels);
            DECLARE_PARAMETER(int,option);

    };

    class Sinusoidal_Volume_Grating : public Grating {
        protected:
            void setup();
        public:
            virtual COMPLEX fourierx(int order,int level,int recip=0);
            virtual COMPLEX eps(double x, int level, int direction);

        private:
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,thick);
            DECLARE_PARAMETER(double,tilt);
            DECLARE_PARAMETER(dielectric_function,minimum);
            DECLARE_PARAMETER(dielectric_function,maximum);
            DECLARE_PARAMETER(int,nlevels);
    };


    class Overlaid_Grating : public Grating {
        public:
            virtual COMPLEX fourierx(int order,int level,int recip=0);
            virtual COMPLEX fouriery(int order,int level,int recip=0);
            virtual COMPLEX fourierz(int order,int level,int recip=0);
            virtual COMPLEX fouriermux(int order,int level,int recip=0);
            virtual COMPLEX fouriermuy(int order,int level,int recip=0);
            virtual COMPLEX fouriermuz(int order,int level,int recip=0);

        protected:
            void setup();

        private:
            DECLARE_MODEL();
            DECLARE_PARAMETER(Grating_Ptr,bottom);
            DECLARE_PARAMETER(Grating_Ptr,top);
            DECLARE_PARAMETER(double,overlay);
            DECLARE_PARAMETER(double,separation);
    };


}
#endif
