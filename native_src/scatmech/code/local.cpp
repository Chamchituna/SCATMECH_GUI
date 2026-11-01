











#include <complex>
#include "scatmech.h"
#include "brdf.h"
#include "local.h"
#include "askuser.h"

using namespace std;


namespace SCATMECH {









    JonesMatrix
    Local_BRDF_Model::
    JonesDSC(double thetai,double thetas,double phis,double rotation, Coordinate_System cs)
    {
        set_geometry(thetai,thetas,phis,rotation);
        JonesMatrix j = jonesDSC();
        convert(j,cs);
        return j;
    }








    MuellerMatrix
    Local_BRDF_Model::
    MuellerDSC(double thetai,double thetas,double phis,double rotation,Coordinate_System cs)
    {
        set_geometry(thetai,thetas,phis,rotation);
        MuellerMatrix m = muellerDSC();
        convert(m,cs);
        return m;
    }

    void
    Local_BRDF_Model::
    setup()
    {
        BRDF_Model::setup();
    }

    JonesMatrix
    Local_BRDF_Model::
    JonesDSC(const Vector& source, const Vector& viewer, const Vector& normal, const Vector& xaxis,
             Coordinate_System cs)
    {
        if (set_geometry(source,viewer,normal,xaxis)) {
            return JonesDSC(thetai,thetas,phis,rotation,cs);
        } else {
            return JonesZero();
        }
    }

    MuellerMatrix
    Local_BRDF_Model::
    MuellerDSC(const Vector& source, const Vector& viewer, const Vector& normal, const Vector& xaxis,
               Coordinate_System cs)
    {
        if (set_geometry(source,viewer,normal,xaxis)) {
            return MuellerDSC(thetai,thetas,phis,rotation,cs);
        } else {
            return MuellerZero();
        }
    }

    DEFINE_VIRTUAL_MODEL(Local_BRDF_Model,BRDF_Model,
                         "All models for discrete defects.");

    DEFINE_PARAMETER(Local_BRDF_Model,double,density,"Density [um^-2]","1",0xFF);


}

