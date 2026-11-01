














#include "filmtran.h"

using namespace std;
using namespace SCATMECH;


int main()
{
    try {


        double lambda = AskUser("Wavelength",0.633);


        dielectric_function substrate =
            dielectric_function::AskUser("Substrate",optical_constant(1.5,0));


        dielectric_stack stack = dielectric_stack::AskUser("Stack","");


        for (double theta=0; theta<90; theta+=1) {


            cout << theta << tab
                 << norm(stack.rp12(theta*deg,lambda,vacuum,substrate)) << tab
                 << norm(stack.rs12(theta*deg,lambda,vacuum,substrate)) << endl;
        }

    } catch (exception& e) {

        cerr << e.what() << endl;

    }

    return 0;
}
