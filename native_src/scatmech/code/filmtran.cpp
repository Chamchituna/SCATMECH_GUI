


































#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include "scatmech.h"
#include "filmtran.h"
#include "fresnel.h"
#include "askuser.h"

using namespace std;


namespace SCATMECH {







    static const COMPLEX I = COMPLEX(0,1);







    Thin_Film_Transfer_Matrix
    Transfer_Matrix(const optical_constant& n,double thickness,COMPLEX theta,double lambda)
    {
        Thin_Film_Transfer_Matrix tftm;
        COMPLEX cos_thetai=cos_internal_angle(n,theta);
        COMPLEX klcostheta = (COMPLEX)n*2.*pi*thickness/lambda*cos_thetai;
        COMPLEX p = cos_thetai*(COMPLEX)n;
        COMPLEX q = cos_thetai/(COMPLEX)n;
        COMPLEX sinklcostheta = sin(klcostheta);
        tftm.As = tftm.Ap = cos(klcostheta);
        tftm.Bs = -I/p*sinklcostheta;
        tftm.Cs = -I*p*sinklcostheta;
        tftm.Bp = -I/q*sinklcostheta;
        tftm.Cp = -I*q*sinklcostheta;
        tftm.Ds = tftm.Dp = tftm.As;
        return tftm;
    }







    COMPLEX
    Thin_Film_Transfer_Matrix::
    rs(COMPLEX theta,const optical_constant& n0,const optical_constant& nt)
    {
        COMPLEX p0 = (COMPLEX)n0*cos_internal_angle(n0,theta);
        COMPLEX pt = (COMPLEX)nt*cos_internal_angle(nt,theta);
        return ((As+Bs*pt)*p0-(Cs+Ds*pt))/((As+Bs*pt)*p0+(Cs+Ds*pt));
    }

    COMPLEX
    Thin_Film_Transfer_Matrix::
    rp(COMPLEX theta,const optical_constant& n0,const optical_constant& nt)
    {
        COMPLEX q0 = cos_internal_angle(n0,theta)/(COMPLEX)n0;
        COMPLEX qt = cos_internal_angle(nt,theta)/(COMPLEX)nt;
        return ((Ap+Bp*qt)*q0-(Cp+Dp*qt))/((Ap+Bp*qt)*q0+(Cp+Dp*qt));
    }

    COMPLEX
    Thin_Film_Transfer_Matrix::
    ts(COMPLEX theta,const optical_constant& n0,const optical_constant& nt)
    {
        COMPLEX p0 = (COMPLEX)n0*cos_internal_angle(n0,theta);
        COMPLEX pt = (COMPLEX)nt*cos_internal_angle(nt,theta);
        return (2.*p0)/((As+Bs*pt)*p0+(Cs+Ds*pt));
    }

    COMPLEX
    Thin_Film_Transfer_Matrix::
    tp(COMPLEX theta,const optical_constant& n0,const optical_constant& nt)
    {
        COMPLEX p0 = cos_internal_angle(n0,theta)/(COMPLEX)n0;
        COMPLEX pt = cos_internal_angle(nt,theta)/(COMPLEX)nt;
        return (2.*p0)/((Ap+Bp*pt)*p0+(Cp+Dp*pt))*(COMPLEX)n0/(COMPLEX)nt;
    }





    Thin_Film_Transfer_Matrix
    Thin_Film_Transfer_Matrix::
    operator*(Thin_Film_Transfer_Matrix a)
    {
        Thin_Film_Transfer_Matrix t;
        t.As = As*a.As+Bs*a.Cs;
        t.Ap = Ap*a.Ap+Bp*a.Cp;
        t.Bs = As*a.Bs+Bs*a.Ds;
        t.Bp = Ap*a.Bp+Bp*a.Dp;
        t.Cs = Cs*a.As+Ds*a.Cs;
        t.Cp = Cp*a.Ap+Dp*a.Cp;
        t.Ds = Cs*a.Bs+Ds*a.Ds;
        t.Dp = Cp*a.Bp+Dp*a.Dp;
        return t;
    }

























    Thin_Film_Transfer_Matrix
    dielectric_stack::
    matrix_forward(COMPLEX theta,double lambda) const
    {
        Thin_Film_Transfer_Matrix M;
        int i;
        for (i=n-1; i>=0; --i)
            M = M*Transfer_Matrix(e[i].index(lambda),t[i],theta,lambda);
        return M;
    }




    Thin_Film_Transfer_Matrix
    dielectric_stack::
    matrix_backward(COMPLEX theta,double lambda) const
    {
        Thin_Film_Transfer_Matrix M;
        int i;
        for (i=0; i<n; ++i)
            M = M*Transfer_Matrix(e[i].index(lambda),t[i],theta,lambda);
        return M;
    }

	void dielectric_stack::reverse()
	{
		for (int i=0;i<n/2;++i) {
			swap(e[i],e[n-i-1]);
			swap(t[i],t[n-i-1]);
		}
	}









    COMPLEX
    dielectric_stack::
    rs12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        Thin_Film_Transfer_Matrix M(matrix_forward(theta,lambda));
        return M.rs(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    rp12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        Thin_Film_Transfer_Matrix M(matrix_forward(theta,lambda));
        return M.rp(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    ts12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        Thin_Film_Transfer_Matrix M(matrix_forward(theta,lambda));
        return M.ts(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    tp12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        Thin_Film_Transfer_Matrix M(matrix_forward(theta,lambda));
        return M.tp(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    rs21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        Thin_Film_Transfer_Matrix M(matrix_backward(theta,lambda));
        return M.rs(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    rp21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        Thin_Film_Transfer_Matrix M(matrix_backward(theta,lambda));
        return M.rp(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    ts21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        Thin_Film_Transfer_Matrix M(matrix_backward(theta,lambda));
        return M.ts(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    tp21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        Thin_Film_Transfer_Matrix M(matrix_backward(theta,lambda));
        return M.tp(theta,n0.index(lambda),nt.index(lambda));
    }

    JonesMatrix
    dielectric_stack::
    r12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        return JonesMatrix(rp12(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           rs12(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           (COMPLEX)(0.),
                           (COMPLEX)(0.));
    }

    JonesMatrix
    dielectric_stack::
    t12(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        return JonesMatrix(tp12(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           ts12(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           (COMPLEX)(0.),
                           (COMPLEX)(0.));
    }

    JonesMatrix
    dielectric_stack::
    r21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        return JonesMatrix(rp21(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           rs21(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           (COMPLEX)(0.),
                           (COMPLEX)(0.));
    }

    JonesMatrix
    dielectric_stack::
    t21(COMPLEX theta,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        return JonesMatrix(tp21(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           ts21(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           (COMPLEX)(0.),
                           (COMPLEX)(0.));
    }

    static
    COMPLEX
    ArcSin(COMPLEX a)
    {
        if (imag(a)==0.) a += COMPLEX(0,-1E-10);
        COMPLEX temp = sqrt(1. - sqr(a)) + COMPLEX(0,1)*a;
        return COMPLEX(arg(temp),-log(abs(temp)));
    }

    COMPLEX
    dielectric_stack::
    rs12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        Thin_Film_Transfer_Matrix M(matrix_forward(theta,lambda));
        return M.rs(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    rp12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        Thin_Film_Transfer_Matrix M(matrix_forward(theta,lambda));
        return M.rp(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    ts12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        Thin_Film_Transfer_Matrix M(matrix_forward(theta,lambda));
        return M.ts(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    tp12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        Thin_Film_Transfer_Matrix M(matrix_forward(theta,lambda));
        return M.tp(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    rs21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        Thin_Film_Transfer_Matrix M(matrix_backward(theta,lambda));
        return M.rs(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    rp21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        Thin_Film_Transfer_Matrix M(matrix_backward(theta,lambda));
        return M.rp(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    ts21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        Thin_Film_Transfer_Matrix M(matrix_backward(theta,lambda));
        return M.ts(theta,n0.index(lambda),nt.index(lambda));
    }

    COMPLEX
    dielectric_stack::
    tp21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        Thin_Film_Transfer_Matrix M(matrix_backward(theta,lambda));
        return M.tp(theta,n0.index(lambda),nt.index(lambda));
    }

    JonesMatrix
    dielectric_stack::
    r12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        return JonesMatrix(rp12(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           rs12(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           (COMPLEX)(0.),
                           (COMPLEX)(0.));
    }

    JonesMatrix
    dielectric_stack::
    t12i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        return JonesMatrix(tp12(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           ts12(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           (COMPLEX)(0.),
                           (COMPLEX)(0.));
    }

    JonesMatrix
    dielectric_stack::
    r21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        return JonesMatrix(rp21(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           rs21(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           (COMPLEX)(0.),
                           (COMPLEX)(0.));
    }

    JonesMatrix
    dielectric_stack::
    t21i(COMPLEX angle,double lambda,const dielectric_function& n0,const dielectric_function& nt) const
    {
        COMPLEX theta = ArcSin(sin(angle)*(COMPLEX)n0.index(lambda));
        return JonesMatrix(tp21(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           ts21(theta,lambda,n0.index(lambda),nt.index(lambda)),
                           (COMPLEX)(0.),
                           (COMPLEX)(0.));
    }













    dielectric_stack
    dielectric_stack::
    AskUser(const string& prompt,const string& deflt)
    {

        dielectric_stack stack;

        while (1) {
            stack.wash();

            try {

                string prompt2 = prompt
                                 + "\n(':filename' or index and thickness) <"
                                 + deflt + "> :";
                SCATMECH_output << prompt2;


                string response = SCATMECH_input.getstrline();


                if (response.size()==0) {

                    stack.read_stack_file(deflt);
                } else if (response[0]==':') {

                    istringstream responsestr(response);
                    responsestr.get();
                    string filename;
                    responsestr >> filename;
                    stack.read_stack_file(filename);
                } else {
                    int i=1;
                    bool finished=false;

                    string dielfunc;
                    double thickness=-1;

                    while (!finished) {


                        istringstream sresponse(response);
                        dielfunc.clear();
                        thickness=-1;
                        sresponse >> dielfunc >> thickness;


                        if (dielfunc.size()==0) finished=true;

                        else {

                            if (thickness<0.) throw SCATMECH_exception("Thickness must be greater than zero.");


                            stack.grow(dielectric_function(dielfunc),thickness);


                            ostringstream message;
                            message << "Layer #" << ++i << ":";
                            SCATMECH_output << message.str() << endl;


                            response = getstrline(sresponse);
                            string next;
                            istringstream temp(response);
                            temp >> next;
                            if (next.size()==0) {
                                response = SCATMECH_input.getstrline();
                            }
                        }
                    }
                }
                ostringstream message;
                stack.print(message);
                SCATMECH_output << message.str();
                return stack;
            }
            catch(const exception& e) {
                SCATMECH_output << e.what() << endl;

            }
        }
    }

    void
    dielectric_stack::
    read_stack_file(const string& film_file_name)
    {

        if (film_file_name.size()==0) {
            wash();
            return;
        }


        string fname = find_file(film_file_name);
        ifstream_with_comments film_file(fname.c_str());


        if (!film_file) {
            throw SCATMECH_exception("Cannot open file: " + film_file_name);
        }

        do {

            string material;
            film_file >> material;


            if (material.size()==0) return;


            double thickness=-1.;
            film_file >> thickness;
            if (thickness<0.) throw SCATMECH_exception("Thickness must be greater than zero.");


            grow(dielectric_function(material),thickness);

        } while (!film_file.eof());
    }




    void
    dielectric_stack::
    print(ostream& os) const
    {
        os << endl << "#  index                 thickness " << endl
           << "-----------------------------------" << endl;
        if (n==0) {
            os << "No layers" << endl;
        } else {
            for (int i=0; i<n; ++i) {
                string newname = e[i].get_name().substr(0,20);
                while (newname.size()<20) newname += ' ';
                os << i << ": "
                   << newname
                   << "  " << t[i] << endl;
            }
        }
    }

    ostream&
    operator<<(std::ostream& os,const dielectric_stack& ds)
    {
        if (ds.n==0) {
            os << "(empty)";
        } else {
            for (int i=0; i<ds.n; ++i) {
                string name = ds.e[i].get_name();
                os << name << ' ' << ds.t[i];
                if (i!=ds.n-1) os << ' ';
            }
        }
        return os;
    }

    template <>
    void
    ModelParameterSet(dielectric_stack& variable,const string& parameter,const string& value)
    {
        if (parameter.empty()) {
            variable.wash();
			istringstream iss(value);
            string material;
            iss >> material;
            if (iss.fail()) return;
            double thickness;
            iss >> thickness;
            if (iss.fail()) {
                variable.read_stack_file(material);
                return;
            }
            variable.grow(material,thickness);
            while (1) {
                iss >> material;
                if (iss.fail()) return;
                iss >> thickness;
                if (iss.fail()) throw SCATMECH_exception("Invalid thickness parameter for dielectric_stack");
                variable.grow(material,thickness);
            }
        } else if (parameter == "grow") {
            istringstream iss(value);
            string material;
            double thickness;
            iss >> material >> thickness;
            if (iss.fail()) throw SCATMECH_exception("Invalid parameter for dielectric_stack");
            variable.grow(material,thickness);
        } else if (parameter == "wash") {
            variable.wash();
        }
    }

    template <>
    string
    ModelParameterGet(dielectric_stack& variable,const string& parameter)
    {
        ostringstream oss;
        oss << variable;
        return oss.str();
    }

    template <>
    void
    ModelParameterAskUser(dielectric_stack& value,const string& prompt)
    {
        value = dielectric_stack::AskUser(prompt);
    }

    void Register(const StackModel* x)
    {
        static bool Models_Registered = false;
        if (!Models_Registered) {
            Models_Registered=true;

            Register_Model(StackModel);
			Register_Model(No_StackModel);
            Register_Model(Stack_StackModel);
            Register_Model(SingleFilm_StackModel);
			Register_Model(DoubleFilm_StackModel);
            Register_Model(GradedFilm_StackModel);

        }
    }

	DEFINE_VIRTUAL_MODEL(StackModel,Model,"Abstract model for film stacks, allowing parametric variation of dielectric stacks");

	DEFINE_MODEL(No_StackModel,StackModel,"An empty film stack");

	DEFINE_MODEL(Stack_StackModel,StackModel,"StackModel that takes a dielectric stack");
	DEFINE_PARAMETER(Stack_StackModel,dielectric_stack,stack,"Film stack","",0xFF);

	DEFINE_MODEL(SingleFilm_StackModel,StackModel,"A single film");
	DEFINE_PARAMETER(SingleFilm_StackModel,dielectric_function,material,"Material","(1.5,0)",0xFF);
	DEFINE_PARAMETER(SingleFilm_StackModel,double,thickness,"Thickness [um]","0.1",0xFF);

	DEFINE_MODEL(DoubleFilm_StackModel,StackModel,"A double film");
	DEFINE_PARAMETER(DoubleFilm_StackModel,dielectric_function,material1,"Material deposited first","(1.5,0)",0xFF);
	DEFINE_PARAMETER(DoubleFilm_StackModel,double,thickness1,"Thickness of first film [um]","0.1",0xFF);
	DEFINE_PARAMETER(DoubleFilm_StackModel,dielectric_function,material2,"Material deposited second","(1.5,0)",0xFF);
	DEFINE_PARAMETER(DoubleFilm_StackModel,double,thickness2,"Thickness of second film[um]","0.1",0xFF);

	DEFINE_MODEL(GradedFilm_StackModel,StackModel,"A graded index film");
	DEFINE_PARAMETER(GradedFilm_StackModel,dielectric_function,start,"Material closest to substrate","(1.5,0)",0xFF);
	DEFINE_PARAMETER(GradedFilm_StackModel,dielectric_function,end,"Material farthest from substrate","(1,0)",0xFF);
	DEFINE_PARAMETER(GradedFilm_StackModel,double,thickness,"Thickness [um]","0.1",0xFF);
	DEFINE_PARAMETER(GradedFilm_StackModel,int,steps,"Number of steps","5",0xFF);

}



