











#ifndef SCATMECH_SCATEVAL_H
#define SCATMECH_SCATEVAL_H

#include "scatmech.h"
#include <string>
#include <stack>
#include <vector>
#include <sstream>
#include <map>

namespace SCATMECH {






    class Evaluator {
        public:


            typedef std::vector<double> VDOUBLE;


            typedef std::string STRING;


            typedef std::map<STRING,double> VMAP;


            Evaluator(const STRING& str,
                      const VMAP& _variables = VMAP(),
                      bool _top=true
                     ) : input(str), variables(_variables), top(_top) {
				if (top) {

					try {
						evaluate();
					}
					catch (std::exception& e) {
						throw SCATMECH_exception("Error during evaluation of expression \"" + str + "\": " + e.what());
					}
				} else {
					evaluate();
				}
            }


            operator double() const {
                if (result.size()!=1) error("Call to double() with more than one value");
                return result[0];
            }


            STRING ResultString() const;


            double Result(int i) {
                return result[i];
            }


            int NResult() {
                return result.size();
            }

        private:

            const VMAP& variables;

            std::istringstream input;
            std::stack<double> val_stack;
            std::stack<int> op_stack;
            std::stack<int> prec_stack;
            VDOUBLE result;

            void evaluate();
            void operate();
            double call_function(const STRING& s,const VDOUBLE& args);
            void error(const std::string& message) const;
            bool lower_prec(int prec);
            void get_value();
            void get_operator();
            VDOUBLE get_paren(bool allow_empty=false);

            bool top;
    };

    void Register_Evaluator_Function(double (*f)(const std::vector<double>& args),const std::string& name,int nargs);
    bool Evaluate_Function(double *value,const std::string& name, const std::vector<double>& args);
    std::vector<double> OneLineFunction(const std::string& expression,double x,Evaluator::VMAP variables=Evaluator::VMAP());
}

#endif
