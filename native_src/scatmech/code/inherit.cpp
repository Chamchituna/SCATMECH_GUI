











#include "scatmech.h"
#include "inherit.h"
#include "askuser.h"
#include <sstream>

#include "filmtran.h"
#include "scattabl.h"
#include "dielfunc.h"

using namespace std;


namespace SCATMECH {


    bool Inheritance::show_model_mode=false;
    int Model::quiet=0;





    void stringlist(const string& s,deque<string>& result)
    {
        result.clear();
        result.push_back("");

        for (string::const_iterator ps=s.begin(); ps!=s.end(); ++ps) {
            if (isspace(*ps)) {
                if (!result.back().empty()) result.push_back("");
            } else {
                result.back() += *ps;
            }
        }
        if (result.back().empty()) result.pop_back();
    }




    void Wrap_String(const string& in,string& out,int indent,int width)
    {
        deque<string> ls;
        stringlist(in,ls);

        string indenter="\n";
        for (int i=0; i<indent; ++i) indenter += space;
        int wide = width-indent;

        int currentwidth = 0;
        for (deque<string>::const_iterator pls=ls.begin(); pls!=ls.end(); ++pls) {
            if (((int)pls->size()+currentwidth+1)>=wide) {
                out += indenter;
                currentwidth = pls->size()+1;
            } else {
                currentwidth += pls->size()+1;
            }
            out += *pls;
            out += space;
        }
    }





    void
    display_descriptions(const Inheritance* inherit,const InheritanceList& ilist)
    {




        SCATMECH_output << endl << "-- Available " << inherit->get_name() << " Classes --" << endl;


        InheritanceList::const_iterator v;



        int jj = 0,i;


        for (v = ilist.begin(); v != ilist.end(); ++v) {
            int k = ((*v)->get_name()).size();
            if ((*v)->children.size()>1 && (*v)!=inherit) {
                k+=7;
            }
            if (k>jj) jj=k;
        }

        for (i=0,v = ilist.begin(); v != ilist.end(); ++i,++v) {
            int k = ((*v)->get_name()).size();
            ostringstream message;
            message << "(" << (char)(i+'A') << ") ";
            SCATMECH_output << message.str();
            int menu=0;
            if ((*v)->children.size()>0 && (*v)!=inherit) {
                SCATMECH_output << "<MENU> ";
                menu=7;
            }
            SCATMECH_output << (*v)->get_name();
            if (Inheritance::show_model_mode) {
                for (int j=0; j<jj-k+2-menu; ++j) SCATMECH_output << '.';
                string message;
                Wrap_String((*v)->get_description(),message,jj+6,80);
                SCATMECH_output << message;
            }
            SCATMECH_output << endl;
        }
    }







    Model*
    Inheritance::
    get_model(bool top) const
    {



        while (1) {

            InheritanceList ilist;

            get_model_inheritance_list(ilist,2);



            if (ilist.size()==1) return ilist[0]->make();

            InheritanceList::iterator v;
            for (v=ilist.begin(); v!=ilist.end(); ++v) {

                if ((*v)->is_virtual()&&((*v)->children.size()==0||(*v)==this)) {
                    v = ilist.erase(v);
                }
            }


            display_descriptions(this,ilist);


            SCATMECH_output << endl;
            if (show_model_mode) {
                SCATMECH_output << "(?) Turn off descriptions" << endl;
            } else {
                SCATMECH_output << "(?) Turn on descriptions" << endl;
            }

            if (!top) SCATMECH_output << "(^) Return to previous menu" << endl;

            string response = AskUser("\nChoose one",string("A"));

            Model* model = get_named_model(response,true);

            if (model) {
                return model;
            } else if (response[0]=='?') {
                show_model_mode= !show_model_mode;
            } else if (response[0]=='^' && !top) {
                return (Model*)0;
            } else if (response.size()==1) {

                int k = toupper(response[0])-'A';

                if (k>=0&&k<(int)ilist.size()) {

                    InheritanceList::const_iterator v = ilist.begin();
                    advance(v,k);
                    const Inheritance* ib = *v;
                    if (ib==this) {
                        return ib->make();
                    } else {
                        Model* model = ib->get_model(false);

                        if (model) {
                            return model;
                        } else {
                            return 0;
                        }
                    }
                }
            }
        }
    }




    void
    Inheritance::
    get_parameters(ModelParameterList& result,bool top) const
    {
        if (top) result.clear();
        if (parent) {
            parent->get_parameters(result,false);
        }
        for (ModelParameterList::const_iterator p=parameters.begin(); p!=parameters.end(); ++p) {
            result.push_back(*p);
        }
    }





    const  ModelParameterBase*
    Inheritance::get_modelparameter(const std::string& param) const
    {
        ModelParameterList mpl = get_parameters();
        for (ModelParameterList::const_iterator p=parameters.begin(); p!=parameters.end(); ++p) {
            if ((*p)->name == param) return (*p);
        }
        return NULL;
    }





    void
    Inheritance::
    get_model_inheritance_list(InheritanceList& result,int option) const
    {
        if (option) {
            result.clear();
            result.push_back(this);
        }

        for (InheritanceList::const_iterator v=children.begin(); v != children.end(); ++v) {
            result.push_back(*v);
            if (option!=2) (*v)->get_model_inheritance_list(result,0);
        }
    }




    Model*
    Inheritance::
    get_named_model(const std::string& model,bool nothrow) const
    {
        const Inheritance *inheritance = get_named_inheritance(model,nothrow);
        if (inheritance) {
            if (!inheritance->virt)    return inheritance->make();
        }
        if (!nothrow) throw SCATMECH_exception("Model " + model + " is not instantiable");
        return (Model*)0;
    }

    const Inheritance*
    Inheritance::
    get_named_inheritance(const std::string& model,bool nothrow) const
    {
		if (model == name) return this;

        InheritanceList ilist;
        get_model_inheritance_list(ilist);

        for (InheritanceList::const_iterator v=ilist.begin(); v!=ilist.end(); ++v) {

            if (model==(*v)->get_name()) return *v;
        }


        if (!nothrow) {
            ostringstream os;
            os << string("Cannot find a " + name + " named " + model) << endl;
            os << "Valid choices are" << endl;
            for (InheritanceList::const_iterator v = ilist.begin(); v!=ilist.end(); ++v) {
                os << "   " << (*v)->get_name() << endl;
            }
            throw SCATMECH_exception(os.str());
        }  else {
            return (const Inheritance* )0;
        }
	}

    ModelParameterBase::
    ModelParameterBase(const std::string& _name,
                       const std::string& _description,
                       const std::string& _type,
                       const std::string& _defaultvalue,
                       int _recalclevel,
                       Inheritance& inherit)
        : ParameterInfo(_name,_description,_type,_defaultvalue,_recalclevel)
    {
        inherit.add_parameter(this);
    }

    void
    Model::
    set_parameter_base(const std::string& parameter,const std::string& value)
    {
        if (parameter=="recalc") {
            recalc |= from_string<int>(value);
        } else if (parameter=="quiet") {
            int _quiet = from_string<int>(value);
            set_quiet(_quiet);
        } else {
            ModelParameterList mplist = get_inheritance().get_parameters();

            string _parameter = parameter;
            string _subparameter="";
            int pos = _parameter.find('.');

            if (pos!=string::npos) {
                _parameter = parameter.substr(0,pos);
                _subparameter = parameter.substr(pos+1,parameter.size());
                if (_subparameter=="")
                    error("Subparameter missing: " + parameter);
            }

            ModelParameterList::iterator q;
            for (q = mplist.begin(); q != mplist.end(); ++q) {
                const ModelParameterBase *param = *q;
                string name = param->name;

                if (_parameter==name) {
                    if (value!="default") {
                        param->set_parameter(this,_subparameter,value);
                    } else {
                        if (param->defaultvalue!="") {
                            param->set_parameter(this,_subparameter,param->defaultvalue);
                        }
                    }
                    return;
                }
            }

            error("Unknown parameter " + parameter + " in set_parameter_base");
        }
    }

    string
    Model::
    get_parameter_base(const string& parameter) const
    {
        ostringstream result;

        if (parameter=="recalc") {
            result << recalc;
        } else if (parameter=="quiet") {
            result << get_quiet();
        } else {
            ModelParameterList mplist = get_inheritance().get_parameters();

            string _parameter = parameter;
            string _subparameter="";
            int pos = _parameter.find('.');

            if (pos!=string::npos) {
                _parameter = parameter.substr(0,pos);
                _subparameter = parameter.substr(pos+1,parameter.size());
                if (_subparameter=="")
                    error("Subparameter missing: " + parameter);
            }

            ModelParameterList::iterator q;
            for (q = mplist.begin(); q != mplist.end(); ++q) {
                const ModelParameterBase *param = *q;
                string name = param->name;

                if (_parameter==name) {
                    result << param->get_parameter(this,_subparameter);
                    return result.str();
                }
            }

            error("Unknown parameter " + parameter + " in get_parameter");
        }
        return result.str();
    }

    ParameterInfo
    Model::
    get_parameter_info(const string& parameter) const
    {
        ParameterInfo result;
        if (parameter=="recalc") {
            result.name="recalc";
            result.type="int";
        } else if (parameter=="quiet") {
            result.name="quiet";
            result.type="int";
        } else {
            ModelParameterList mplist = get_inheritance().get_parameters();

            string _parameter = parameter;
            string _subparameter="";
            int pos = _parameter.find('.');

            if (pos!=string::npos) {
                _parameter = parameter.substr(0,pos);
                _subparameter = parameter.substr(pos+1,parameter.size());
                if (_subparameter=="")
                    error("Subparameter missing: " + parameter);
            }

            ModelParameterList::iterator q;
            for (q = mplist.begin(); q != mplist.end(); ++q) {
                const ModelParameterBase *param = *q;
                string name = param->name;

                if (_parameter==name) {
                    if (_subparameter.size()==0) {
                        result = *param;
                    } else {
                        if (param->get_inheritance(NULL)!=NULL) {
                            Model* m =(Model*)(param->get_ptr(this));
                            result = m->get_parameter_info(_subparameter);
                        } else {
                            error("Parameter " + parameter + " does not have subparameters");
                        }
                    }
                    return result;
                }
            }

            error("Unknown parameter " + parameter + " in get_parameter_info");
        }
        return result;
    }

    void
    Model::
    init()
    {
        ModelParameterList mplist = get_inheritance().get_parameters();

        ModelParameterList::iterator q;
        for (q = mplist.begin(); q != mplist.end(); ++q) {
            const ModelParameterBase *param = *q;
            if (param->defaultvalue!="") {
                param->set_parameter(this,"",param->defaultvalue);
            }
        }
    }

    void
    Model::
    AskUser()
    {
        init();

        ModelParameterList mplist = get_inheritance().get_parameters();

        ModelParameterList::iterator q;
        for (q = mplist.begin(); q != mplist.end(); ++q) {
            const ModelParameterBase* param = *q;
            param->AskUser(this);
        }
    }

    void
    Model::
    get_parameter_names_base(StringList& plist, const string& prefix) const
    {
        ModelParameterList mpl = get_inheritance().get_parameters();

        for (ModelParameterList::iterator q=mpl.begin(); q!=mpl.end(); ++q) {
            plist.push_back(prefix+(*q)->name);
            if ((*q)->get_inheritance(this)!=NULL) {
                Model* m = (Model*)((*q)->get_ptr(this));
                m->get_parameter_names_base(plist,prefix+(*q)->name+'.');
            }
        }
    }

    void
    Model::
    print_parameters_base(ostream& os,const string& prefix) const
    {









        ModelParameterList mplist = get_inheritance().get_parameters();
        for (ModelParameterList::iterator p = mplist.begin(); p!=mplist.end(); ++p) {
            const ModelParameterBase *param = *p;
            string name = param->name;
            string type = param->type;
            string description = param->description;
			string newprefix = (prefix != "") ? (prefix + ".") : ("");

            if (param->get_inheritance(this)==NULL) {
                os << newprefix << name << " = " << param->get_parameter(this,"")
                   << " (" << type << ": "
                   << description << ")" << endl;
            } else {
                os << newprefix << name << " = " << param->get_parameter(this,"")
                   << " (" << type << ": "
                   << description << ")" << endl;

                ((Model*)(param->get_ptr(this)))->print_parameters_base(os,newprefix+name);
            }
        }
    }

    void
    Model::
    error(const string& message) const
    {
        string message2 = "Error in class " + get_inheritance().get_name() + ": " + message + '\n';
        ostringstream out;
        print_parameters(out);
        throw SCATMECH_exception(message2+out.str());
    }

    Inheritance Model::inheritance("Model","Generalized Model","",0,0);


}


