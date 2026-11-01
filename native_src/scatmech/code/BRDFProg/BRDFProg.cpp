













#include "brdf.h"
#include <fstream>

using namespace std;
using namespace SCATMECH;

int main(int argc, char** argv)
{
    try {

        istream* in = &cin;
        if (argc == 3 && string(argv[1]) == "--batch") {
            in = new ifstream(argv[2]);
            if (!static_cast<ifstream*>(in)->is_open()) {
                cerr << "Could not open batch file: " << argv[2] << endl;
                return 1;
            }
        }


        streambuf* orig_cin = cin.rdbuf();
        if (in != &cin) cin.rdbuf(in->rdbuf());

        double thetai = AskUser("Incident Angle ", 45.) * deg;
        double thetas_1 = AskUser("Initial Scattering Angle ", 45.) * deg;
        double thetas_2 = AskUser("Final Scattering Angle ", 45.) * deg;
        double thetas_3 = AskUser("Step Scattering Angle ", 1.) * deg;
        double phis_1 = AskUser("Initial Azimuthal Angle ", 0.) * deg;
        double phis_2 = AskUser("Final Azimuthal Angle ", 180.) * deg;
        double phis_3 = AskUser("Step Azimuthal Angle ", 2.) * deg;

        BRDF_Model_Ptr model = Get_Model_Ptr();
        model->AskUser();

        for (double thetas = thetas_1; thetas <= thetas_2; thetas += thetas_3) {
            for (double phis = phis_1; phis <= phis_2; phis += phis_3) {
                MuellerMatrix m = model->Mueller(thetai, thetas, phis, 0.);
                StokesVector s = m * StokesVectorUnitP();
                cout << thetas / deg << tab
                     << phis / deg << tab
                     << s.eta() / deg << tab
                     << s.DOLP() << tab
                     << s.DOCP() << tab
                     << s.DOP() << tab
                     << s.I() << endl;
            }
        }


        if (in != &cin) {
            cin.rdbuf(orig_cin);
            delete in;
        }

    } catch (exception& e) {
        cerr << e.what() << endl;
    }

    return 0;
}
