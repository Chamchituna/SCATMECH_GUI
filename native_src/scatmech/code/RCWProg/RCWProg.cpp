














#include "rcw.h"

using namespace std;
using namespace SCATMECH;

int main()
{
    try {


        RCW_Model RCW;


        RCW.AskUser();


        RCW.set_quiet(1);


        cerr << endl << "Available parameters:" << endl << endl;
        RCW.print_parameters(cerr);
        cerr << endl;


        string parameter=AskUser("Parameter to vary",string("lambda"));
        double begin=AskUser("Initial value",0.200);
        double end=AskUser("Final value",0.800);
        double step=AskUser("Step value",0.005);


        cout << parameter << tab << "alpha" << tab << "beta" << endl;


        for (double value = begin; value <= end; value += step) {


            RCW.set_parameter(parameter,value);


            MuellerMatrix M = RCW.GetIntensity(0);


            StokesVector S = M*StokesVectorUnitLinear(45.*deg);


            cout << value << tab << S[1]/S[0] << tab << S[2]/S[0] << endl;
        }

    } catch (exception& e) {

        cerr << e.what() << endl;

    }
    return 0;
}
