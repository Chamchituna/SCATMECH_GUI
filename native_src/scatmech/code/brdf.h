












#ifndef SCATMECH_BRDF_H
#define SCATMECH_BRDF_H

#include <limits>
#include "scatmech.h"
#include "mueller.h"
#include "dielfunc.h"
#include "inherit.h"
#include "vector3d.h"

namespace SCATMECH {




    class BRDF_Model : public Model {
        public:



            enum Coordinate_System {
                psps,
                xyxy,
                plane,
                undefined
            };



            enum Type {


                Type_DOWNUP = 0,
                Type_DOWNDOWN = 1,
                Type_UPDOWN = 2,
                Type_UPUP = 3
            };



            BRDF_Model() {
                model_cs = psps;
                last_thetai = last_thetas = last_phis = last_rotation = std::numeric_limits<double>::infinity();
                recalc_on_incident_change = recalc_on_scatter_change = 0x00;
            }






            JonesMatrix Jones(double thetai,
                              double thetas,
                              double phis,
                              double rotation,
                              Coordinate_System cs = psps);



            MuellerMatrix Mueller(double thetai,
                                  double thetas,
                                  double phis,
                                  double rotation,
                                  Coordinate_System cs = psps);




            JonesMatrix Jones(const Vector& source,
                              const Vector& viewer,
                              const Vector& normal,
                              const Vector& xaxis=Vector(0,0,0),
                              Coordinate_System cs=plane);



            MuellerMatrix Mueller(const Vector& source,
                                  const Vector& viewer,
                                  const Vector& normal,
                                  const Vector& xaxis=Vector(0,0,0),
                                  Coordinate_System cs=plane);



            Coordinate_System get_model_cs() const {
                return model_cs;
            }


            DECLARE_MODEL();




            DECLARE_PARAMETER(double,lambda);




            DECLARE_PARAMETER(dielectric_function,substrate);









            DECLARE_PARAMETER(int,type);

        protected:


            virtual void setup();





            Coordinate_System model_cs;



            int set_geometry(double _thetai,double _thetas,double _phis, double _rotation);


            int set_geometry(const Vector& source,const Vector& viewer,const Vector& normal,const Vector& xaxis);


            void convert(JonesMatrix& J,Coordinate_System to) const;

            void convert(MuellerMatrix& M,Coordinate_System to) const;


            double thetai;

            double thetas;

            double phis;

            double rotation;







            virtual JonesMatrix jones() {

                error("Attempt to call jones() for depolarizing model.");
                return JonesMatrix();
            }




            virtual MuellerMatrix mueller()
            {

                return this->jones();
            }






            void set_recalc_on_direction_change(int incident_change,int scatter_change) {
                recalc_on_incident_change = incident_change;
                recalc_on_scatter_change = scatter_change;
            };

        public:

            bool is_reflection() const {
                return type==Type_DOWNUP||type==Type_UPDOWN;
            }

            bool is_transmission() const  {
                return type==Type_DOWNDOWN||type==Type_UPUP;
            }

            bool is_forward() const {
                return type==Type_DOWNDOWN||type==Type_DOWNUP;
            }

            bool is_backward() const  {
                return type==Type_UPDOWN||type==Type_UPUP;
            }


            bool is_down_to_up() const  {
                return type==Type_DOWNUP;
            }

            bool is_down_to_down() const  {
                return type==Type_DOWNDOWN;
            }

            bool is_up_to_down() const  {
                return type==Type_UPDOWN;
            }

            bool is_up_to_up() const  {
                return type==Type_UPUP;
            }


            bool is_inup() const  {
                return is_backward();
            }

            bool is_indown() const  {
                return is_forward();
            }

            bool is_outup() const  {
                return type==Type_DOWNUP||type==Type_UPUP;
            }

            bool is_outdown() {
                return type==Type_DOWNDOWN||type==Type_UPDOWN;
            }

        protected:

            void throw_forward() {
                if (is_forward()) error("No forward (type==0 or 1) code for model");
            }

            void throw_backward() {
                if (is_backward()) error("No backward (type==2 or 3) code for model");
            }

            void throw_transmission() {
                if (is_transmission()) error("No transmission (type==1 or 3) code for model");
            }

            void throw_reflection() {
                if (is_reflection()) error("No reflection (type==0 or 2) code for model");
            }

        private:
            double last_thetai;
            double last_thetas;
            double last_phis;
            double last_rotation;
            int recalc_on_incident_change;
            int recalc_on_scatter_change;
    };




    typedef Model_Ptr<BRDF_Model> BRDF_Model_Ptr;

    void Register(const BRDF_Model* x);


}


#endif

