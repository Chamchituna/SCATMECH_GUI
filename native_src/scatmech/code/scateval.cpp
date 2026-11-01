











#include "scateval.h"
#include "scattabl.h"
#include <map>

using namespace std;

namespace SCATMECH {

    namespace {


        static const char binops[] = ",&|><=+-*/^";


        static const int precs[] = {0,1,1,2,2,2,3,3,5,5,6};



        int member(char c,const char* p) {
            for (int i=0; *p; p++,i++) if (*p == c) return i;
            return -1;
        }


        bool isvstart(char c)
        {
            if (isalpha(c)) return true;
            if (c=='#') return true;
            if (c=='$') return true;
            if (c=='@') return true;
            return false;
        }



        bool isvcont(char c)
        {
            if (isvstart(c)) return true;
            if (isalnum(c)) return true;
            if (c=='.') return true;
            if (c=='_') return true;
            return false;
        }

        typedef map< string, double (*)(const std::vector<double>& args) > Evaluator_Function_Map0;
        typedef map< int, Evaluator_Function_Map0 > Evaluator_Function_Map;
        Evaluator_Function_Map Evaluator_Functions;
    }



    vector<double>
    Evaluator::
    get_paren(bool allow_empty)
    {
        if (input.peek()!='(') {
            error("Expected left parentheses: " + input.str());
        }


        input.ignore();


        input >> ws;


        if (allow_empty && input.peek()==')') {
            input.ignore();
            return vector<double>();
        }


        string contents;
        int level=1;
        while (level>0) {
            if (input.eof()) error("Mismatched parentheses");
            char next = input.get();
            if (next==')') --level;
            if (next=='(') ++level;

            if (level>0) contents += next;
        }

        return Evaluator(contents,variables,false).result;
    }


    void
    Evaluator::
    get_value()
    {

        input >> ws;


        char pre = input.peek();
        int sign = 1;
        if (pre == '-') {
			if (val_stack.size()==0) {
				val_stack.push(0.);
				get_operator();
			} else {
				input.get();
				sign = -1;
			}
        } else if (pre == '+') {
            input.get();
        }


        double x;
        input >> x;
        if (!input.fail()) {

            val_stack.push(sign*x);
        } else {

            input.clear();


            if (input.peek()=='(') {

                vector<double> temp = get_paren();

                if (temp.size()<1) error("Empty parentheses not allowed");
                if (temp.size()>1 && !top) error("Embedded parentheses group");
                val_stack.push(sign*temp[0]);
                for (int i=1; i<(int)temp.size(); ++i) result.push_back(sign*temp[i]);
            } else {

                if (input.peek()==')') error("Mismatched parentheses");

                input >> ws;


                if (isvstart(input.peek())) {
                    string a;
                    while (isvcont(input.peek())) {
                        a += input.get();
                    }

                    input >> ws;


                    vector<double> args;
                    if (input.peek()=='(') args = get_paren(true);


                    double vv = call_function(a,args);
                    val_stack.push(sign*vv);

                } else {
                    error("Invalid value");
                }
            }
        }

        input >> ws;
    }

    void
    Evaluator::
    get_operator()
    {

        input >> ws;

        char next = input.peek();
        int prec=member(next,binops);

        if (prec>=0) {
            while (lower_prec(precs[prec])) operate();
            op_stack.push(input.get());
            prec_stack.push(precs[prec]);
        } else error("Undefined operator: " + next);
        input >> ws;
    }

    void
    Evaluator::
    evaluate()
    {
        get_value();

        while (!input.eof()) {
            get_operator();
            get_value();
        }

        while (!op_stack.empty()) {
            operate();
        }

        if (!val_stack.empty()) result.insert(result.begin(),val_stack.top());
    }

    void
    Evaluator::
    operate()
    {
        int op = op_stack.top();
        op_stack.pop();
        prec_stack.pop();

        double y = val_stack.top();
        val_stack.pop();
        double x = val_stack.top();
        val_stack.pop();

        switch (op) {
            case '+':
                val_stack.push(x+y);
                break;
            case '-':
                val_stack.push(x-y);
                break;
            case '*':
                val_stack.push(x*y);
                break;
            case '/':
                val_stack.push(x/y);
                break;
            case '&':
                val_stack.push(x&&y);
                break;
            case '|':
                val_stack.push(x||y);
                break;
            case '<':
                val_stack.push(x<y);
                break;
            case '>':
                val_stack.push(x>y);
                break;
            case '=':
                val_stack.push(x==y);
                break;
            case '^':
                val_stack.push(pow(x,y));
                break;
            case ',':
                val_stack.push(x);
                result.push_back(y);
                break;

            default:
                error("Invalid binary operator");
        }
    }

    double
    Evaluator::
    call_function(const string& s,const vector<double>& args)
    {
        if (s[0]=='@') {
            string filename = s.substr(1);
            int asize  = args.size();
            if (asize<1 || asize>2) error("File function " + s + " requires 1 or 2 arguments");
            int column = asize==2 ? (int)(args[1]) : 2;
            if (column <= 1) error("File function " + s + " requires column number greater than 1");
            Table table(filename,column);
            return table.value(args[0]);
        }

        if (args.size()==0) {
            VMAP::const_iterator it = variables.find(s);
            if (it != variables.end()) return it->second;
        }

        double result;
        if (Evaluate_Function(&result,s,args)) return result;

        Evaluator_Function_Map::iterator q = Evaluator_Functions.find(args.size());
        if (q!=Evaluator_Functions.end()) {
            Evaluator_Function_Map0::iterator p = q->second.find(s);
            if (p!=q->second.end()) {
                return p->second(args);
            }
        }

        error("Invalid function or value: " + s + " having " + to_string(args.size()) + " arguments.");

        return 0;
    }

    bool
    Evaluator::
    lower_prec(int prec)
    {
        if (prec_stack.empty()) return false;
        if (prec<=prec_stack.top()) return true;
        return false;
    }

    string
    Evaluator::
    ResultString() const
    {
        ostringstream out;
        out.precision(16);

        if (result.size()==1) {
            out << result[0];
            return out.str();
        } else {
            out << '(';
            for (int i=0; i<(int)result.size(); ++i) {
                out << result[i];
                if (i!=result.size()-1) out << ',';
            }
            out << ')';
            return out.str();
        }
    }

    void
    Evaluator::
    error(const string& message) const {
        ostringstream m;
        m << "Evaluator: " << message << endl
          << "Expression: \"" << input.str() << "\"";
        throw SCATMECH_exception(m.str());
    }


    void Register_Evaluator_Function(double (*f)(const std::vector<double>& args),const std::string& name,int nargs)
    {
        Evaluator_Functions[nargs][name] = f;
    }

    bool Evaluate_Function(double *value,const std::string& name, const std::vector<double>& args)
    {
        Evaluator_Function_Map::iterator q = Evaluator_Functions.find(args.size());
        if (q!=Evaluator_Functions.end()) {
            Evaluator_Function_Map0::iterator p = q->second.find(name);
            if (p!=q->second.end()) {
                *value = p->second(args);
                return true;
            }
        }
        return false;
    }

    namespace {

        double exp_EF(const std::vector<double>& args) {
            return exp(args[0]);
        }
        double sin_EF(const std::vector<double>& args) {
            return sin(args[0]);
        }
        double cos_EF(const std::vector<double>& args) {
            return cos(args[0]);
        }
        double tan_EF(const std::vector<double>& args) {
            return tan(args[0]);
        }
        double sind_EF(const std::vector<double>& args) {
            return sin(args[0]*deg);
        }
        double cosd_EF(const std::vector<double>& args) {
            return cos(args[0]*deg);
        }
        double tand_EF(const std::vector<double>& args) {
            return tan(args[0]*deg);
        }
        double asin_EF(const std::vector<double>& args) {
            return asin(args[0]);
        }
        double acos_EF(const std::vector<double>& args) {
            return acos(args[0]);
        }
        double atan_EF(const std::vector<double>& args) {
            return atan(args[0]);
        }
        double asind_EF(const std::vector<double>& args) {
            return asin(args[0])/deg;
        }
        double acosd_EF(const std::vector<double>& args) {
            return acos(args[0])/deg;
        }
        double atand_EF(const std::vector<double>& args) {
            return atan(args[0])/deg;
        }
        double sinh_EF(const std::vector<double>& args) {
            return sinh(args[0]);
        }
        double cosh_EF(const std::vector<double>& args) {
            return cosh(args[0]);
        }
        double tanh_EF(const std::vector<double>& args) {
            return tanh(args[0]);
        }
        double log_EF(const std::vector<double>& args) {
            return log(args[0]);
        }
        double log10_EF(const std::vector<double>& args) {
            return log10(args[0]);
        }
        double sqrt_EF(const std::vector<double>& args) {
            return sqrt(args[0]);
        }
        double abs_EF(const std::vector<double>& args) {
            return fabs(args[0]);
        }
        double not_EF(const std::vector<double>& args) {
            return args[0]==0. ?  1. : 0.;
        }

        double atan2_EF(const std::vector<double>& args) {
            return atan2(args[0],args[1]);
        }
        double atan2d_EF(const std::vector<double>& args) {
            return atan2(args[0],args[1])/deg;
        }
        double min_EF(const std::vector<double>& args) {
            return (args[0]<args[1]) ? args[0] : args[1];
        }
        double max_EF(const std::vector<double>& args) {
            return (args[0]>args[1]) ? args[0] : args[1];
        }
        double or_EF(const std::vector<double>& args) {
            return args[0]||args[1];
        }
        double nor_EF(const std::vector<double>& args) {
            return !(args[0]||args[1]);
        }
        double and_EF(const std::vector<double>& args) {
            return args[0]&&args[1];
        }
        double nand_EF(const std::vector<double>& args) {
            return !(args[0]&&args[1]);
        }
        double resqrt_EF(const std::vector<double>& args) {
            return real(sqrt(COMPLEX(args[0],args[1])));
        }
        double imsqrt_EF(const std::vector<double>& args) {
            return imag(sqrt(COMPLEX(args[0],args[1])));
        }
        double resqr_EF(const std::vector<double>& args) {
            return real(sqr(COMPLEX(args[0],args[1])));
        }
        double imsqr_EF(const std::vector<double>& args) {
            return imag(sqr(COMPLEX(args[0],args[1])));
        }
        double gt_EF(const std::vector<double>& args) {
            return args[0]>args[1] ? 1. : 0.;
        }
        double lt_EF(const std::vector<double>& args) {
            return args[0]<args[1] ? 1. : 0.;
        }
        double ge_EF(const std::vector<double>& args) {
            return args[0]>=args[1] ? 1. : 0.;
        }
        double le_EF(const std::vector<double>& args) {
            return args[0]<=args[1] ? 1. : 0.;
        }
        double eq_EF(const std::vector<double>& args) {
            return args[0]==args[1] ? 1. : 0.;
        }
        double ne_EF(const std::vector<double>& args) {
            return args[0]!=args[1] ? 1. : 0.;
        }

        double if_EF(const std::vector<double>& args) {
            return args[0]!=0 ? args[1] : args[2];
        }

        class Register_Evaluator_Functions {
            public:
                Register_Evaluator_Functions() {
                    Register_Evaluator_Function(&exp_EF,"exp",1);
                    Register_Evaluator_Function(&sin_EF,"sin",1);
                    Register_Evaluator_Function(&cos_EF,"cos",1);
                    Register_Evaluator_Function(&tan_EF,"tan",1);
                    Register_Evaluator_Function(&sind_EF,"sind",1);
                    Register_Evaluator_Function(&cosd_EF,"cosd",1);
                    Register_Evaluator_Function(&tand_EF,"tand",1);
                    Register_Evaluator_Function(&asin_EF,"asin",1);
                    Register_Evaluator_Function(&acos_EF,"acos",1);
                    Register_Evaluator_Function(&atan_EF,"atan",1);
                    Register_Evaluator_Function(&asind_EF,"asind",1);
                    Register_Evaluator_Function(&acosd_EF,"acosd",1);
                    Register_Evaluator_Function(&atand_EF,"atand",1);
                    Register_Evaluator_Function(&sinh_EF,"sinh",1);
                    Register_Evaluator_Function(&cosh_EF,"cosh",1);
                    Register_Evaluator_Function(&tanh_EF,"tanh",1);
                    Register_Evaluator_Function(&log_EF,"log",1);
                    Register_Evaluator_Function(&log10_EF,"log10",1);
                    Register_Evaluator_Function(&sqrt_EF,"sqrt",1);
                    Register_Evaluator_Function(&abs_EF,"abs",1);
                    Register_Evaluator_Function(&not_EF,"not",1);

                    Register_Evaluator_Function(&atan2_EF,"atan2",2);
                    Register_Evaluator_Function(&atan2d_EF,"atan2d",2);
                    Register_Evaluator_Function(&min_EF,"min",2);
                    Register_Evaluator_Function(&max_EF,"max",2);
                    Register_Evaluator_Function(&or_EF,"or",2);
                    Register_Evaluator_Function(&nor_EF,"nor",2);
                    Register_Evaluator_Function(&and_EF,"and",2);
                    Register_Evaluator_Function(&nand_EF,"nand",2);
                    Register_Evaluator_Function(&resqrt_EF,"resqrt",2);
                    Register_Evaluator_Function(&imsqrt_EF,"imsqrt",2);
                    Register_Evaluator_Function(&resqr_EF,"resqr",2);
                    Register_Evaluator_Function(&imsqr_EF,"imsqr",2);
                    Register_Evaluator_Function(&gt_EF,"gt",2);
                    Register_Evaluator_Function(&lt_EF,"lt",2);
                    Register_Evaluator_Function(&ge_EF,"ge",2);
                    Register_Evaluator_Function(&le_EF,"le",2);
                    Register_Evaluator_Function(&eq_EF,"eq",2);
                    Register_Evaluator_Function(&ne_EF,"ne",2);

                    Register_Evaluator_Function(&if_EF,"if",3);
                }
        };
        Register_Evaluator_Functions register_evaluation_functions;
    }

    vector<double> OneLineFunction(const string& expression,double x,Evaluator::VMAP variables)
    {
        if (expression.size()==0) return vector<double>();
        int i=0;
        for (i=0; i<(int)expression.size(); ++i) {
            if (!iswspace(expression[i])) break;
        }
        if (expression[i]!='@') return vector<double>();
        ++i;
        string variable;
        for (; i<(int)expression.size(); ++i) {
            if (expression[i]==':') break;
            if (!iswspace(expression[i])) {
                if (expression[i]!=':') variable.push_back(expression[i]);
            }
        }
        if (i==expression.size()) return vector<double>();
        if (i+1==expression.size()) return vector<double>();
        string formula = expression.substr(i+1);
        variables[variable] = x;
        Evaluator eval(formula,variables);
        vector<double> result(eval.NResult());
        for (i=0; i<eval.NResult(); ++i) {
            result[i]=eval.Result(i);
        }
        return result;
    }

}

