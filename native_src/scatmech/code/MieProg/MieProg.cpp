














#include "miescat.h"

using namespace std;
using namespace SCATMECH;

int main(int argv, char** argc)
{
    try {


        double dtheta = AskUser("Step angle (theta, deg)", 1.0) * deg;


        double dphi   = AskUser("Step azimuth (phi, deg)", 10.0) * deg;


        MieScatterer mie;
        mie.AskUser();


        cout << "Qsca = "  << mie.Qsca()
             << "  Qext = " << mie.Qext()
             << "  Qback = " << mie.Qback() << endl;


        cout << endl
             << "Theta" << tab
             << "Phi"   << tab
             << "S11"   << tab
             << "Pol"   << tab
             << "S33"   << tab
             << "S34"   << endl;


        double m00 = MuellerMatrix(mie.s(0.0))[0][0];


        for (double theta = 0.0; theta <= 180.0 * deg + 1e-12; theta += dtheta) {



            MuellerMatrix m = mie.s(theta);


            for (double phi = 0.0; phi < 360.0 * deg - 1e-12; phi += dphi) {
                cout << theta/deg << tab
                     << phi/deg   << tab
                     << m[0][0]/m00 << tab
                     << m[0][1]/m[0][0] << tab
                     << m[2][2]/m[0][0] << tab
                     << m[2][3]/m[0][0]
                     << endl;
            }


            if (fmod(360.0*deg, dphi) > 1e-12) {
                double phi = 360.0 * deg;
                cout << theta/deg << tab
                     << phi/deg   << tab
                     << m[0][0]/m00 << tab
                     << m[0][1]/m[0][0] << tab
                     << m[2][2]/m[0][0] << tab
                     << m[2][3]/m[0][0]
                     << endl;
            }
        }

    } catch (exception& e) {
        cerr << e.what() << endl;
    }

    return 0;
}
