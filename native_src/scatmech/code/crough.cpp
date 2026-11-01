












#include "scatmech.h"
#include "crough.h"

using namespace std;


namespace SCATMECH {

    void
    Correlated_Roughness_BRDF_Model::
    setup()
    {
        Roughness_BRDF_Model::setup();

        model.set_lambda(lambda);
        model.set_type(type);
        model.set_substrate(substrate);
        model.set_psd(psd);
        SingleFilm_StackModel stack;
		stack.set_material(film);
		stack.set_thickness(thickness);
        model.set_stack(stack.clone());
    }



    JonesMatrix
    Correlated_Roughness_BRDF_Model::
    jones()
    {
        SETUP();

        return model.Jones(thetai,thetas,phis,rotation);
    }

    DEFINE_MODEL(Correlated_Roughness_BRDF_Model,Roughness_BRDF_Model,
                 "Scattering from a single dielectric layer on a substrate with equal and correlated roughness.");

    DEFINE_PARAMETER(Correlated_Roughness_BRDF_Model,double,thickness,"Film thickness [um]","0.05",0xFF);

    DEFINE_PARAMETER(Correlated_Roughness_BRDF_Model,dielectric_function,film,"Overlayer","(1.46,0)",0xFF);


}



