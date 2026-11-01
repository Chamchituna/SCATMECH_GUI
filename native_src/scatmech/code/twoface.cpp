











#include "scatmech.h"
#include "twoface.h"

using namespace std;


namespace SCATMECH {

    void
    Two_Face_BRDF_Model::
    setup()
    {
        Roughness_BRDF_Model::setup();

        model.set_lambda(lambda);
        model.set_type(type);
        model.set_substrate(substrate);
        model.set_psd(psd);
		SingleFilm_StackModel stack0;
		stack0.set_material(film);
		stack0.set_thickness(thickness);
        model.set_stack(stack0.clone());
        model.set_this_layer(face-1);
    }




    JonesMatrix
    Two_Face_BRDF_Model::
    jones()
    {
        SETUP();

        return model.Jones(thetai,thetas,phis,rotation);
    }

    DEFINE_MODEL(Two_Face_BRDF_Model,Roughness_BRDF_Model,
                 "Scattering by a single rough interface of a single dielectric film.");


    DEFINE_PARAMETER(Two_Face_BRDF_Model,dielectric_function,film,"Layer","(1.46,0)",0xFF);

    DEFINE_PARAMETER(Two_Face_BRDF_Model,double,thickness,"Thickness [um]","0.05",0xFF);

    DEFINE_PARAMETER(Two_Face_BRDF_Model,int,face,"Interface number (buried:1, exposed:2)","1",0xFF);




}

