











#ifndef SCATMECH_INHERIT_H
#define SCATMECH_INHERIT_H

#include <deque>
#include "scatmech.h"
#include <memory>

#define SMERROR(message)    error(message,__FILE__,__LINE__);

namespace SCATMECH {




    class Inheritance;
    class Model;
    template <class _Model> class Model_Ptr;
    template <class _Model> class Model_Maker;
    class Model_Maker_Base;
    struct ParameterInfo;
    class ModelParameterBase;
    template <class MODEL,class TYPE> class ModelParameter;
    template <class MODEL,class TYPE> class ModelPtrParameter;
    class Get_Model_Ptr;

    template <class _Model> _Model* _Get_Named_Model(const STRING& name);
    template <class _Model> _Model* _Get_Model();

    typedef std::deque<const Inheritance*> InheritanceList;
    typedef std::unique_ptr<ModelParameterBase> ModelParameterBasePtr;
    typedef std::deque<const ModelParameterBase*> ModelParameterList;
    typedef std::deque<STRING> StringList;




    struct ParameterInfo {

        ParameterInfo() {}

        ParameterInfo(
            const STRING& _name,
            const STRING& _description,
            const STRING& _type,
            const STRING& _defaultvalue,
            int _recalclevel
        ) : name(_name), description(_description),
            type(_type), defaultvalue(_defaultvalue),
            recalclevel(_recalclevel) {}

        STRING name;
        STRING description;
        STRING type;
        STRING defaultvalue;
        int recalclevel;

    };






    class ModelParameterBase : public ParameterInfo {
        public:

            ModelParameterBase(
                const STRING& _name,
                const STRING& _description,
                const STRING& _type,
                const STRING& _defaultvalue,
                int _recalclevel,
                Inheritance& inherit

            );


            virtual ~ModelParameterBase() {}


            virtual void set_parameter(
                Model* model,
                const STRING& subparameter,
                const STRING& value
            ) const = 0;



            virtual void* get_ptr(
                const Model* model
            ) const = 0;


            virtual STRING get_parameter(
                const Model* model,
                const STRING& subparameter
            ) const=0;


            virtual void AskUser(
                Model* model
            ) const = 0;




            virtual const Inheritance* get_inheritance(
                const Model*
            ) const {
                return (Inheritance*)0;
            }
    };




    class Model
    {
        public:

            Model() {
                recalc=-1;
                init();
            }


            virtual ~Model() {}




            void init();




            void AskUser();


            static Inheritance inheritance;




            virtual const Inheritance& get_inheritance() const {
                return inheritance;
            }


            void get_parameter_names(
                StringList& plist
            ) const {
                get_parameter_names_base(plist,"");
            }


            ParameterInfo get_parameter_info(
                const STRING& parameter
            ) const;


            template <class SETTYPE>
            void set_parameter(
                const STRING& parameter,
                const SETTYPE& value
            ) {
                set_parameter_base(parameter,to_string(value,16));
            }


            void set_parameter(
                const STRING& parameter,
                const char* value
            ) {
                set_parameter_base(parameter,std::string(value));
            }


            STRING get_parameter(
                const STRING& parameter
            ) const {
                return get_parameter_base(parameter);
            }



            void print_parameters(
                std::ostream& os,
                const std::string& prefix = ""
            ) const {
                print_parameters_base(os,prefix);
            }


            void set_recalc(int _recalc=0xFF) {
                recalc |= _recalc;
            }



            virtual int get_recalc() const {
                return recalc;
            }


            static void set_quiet(int _quiet) {
                static mutex_t mutex;
                mutex.lock();
                quiet=_quiet;
                mutex.unlock();
            }


            static int get_quiet() {
                return quiet;
            }



            void SETUP() {
                if (get_recalc()) setup();
            }

        protected:




            virtual void setup() {
                RECALC = recalc;
                recalc=0;
            }


            void message(
                const STRING& s
            ) const {
                if (!quiet) {
                    SCATMECH_output << s << ret;
                }
            }


            void error(const STRING& message) const;





            virtual void set_parameter_base(
                const STRING& parameter,
                const STRING& value
            );





            virtual STRING get_parameter_base(
                const STRING& parameter
            ) const;






            virtual void get_parameter_names_base(
                StringList& plist,
                const STRING& prefix
            ) const;





            virtual void print_parameters_base(
                std::ostream& os,
                const std::string& prefix = ""
            ) const;





            int RECALC;

        private:



            int recalc;




            static int quiet;

    };










    class Inheritance {
        public:

            Inheritance(
                const STRING& _name,
                const STRING& _desc,
				const STRING& _parentname,
                Model_Maker_Base *_maker,
                Inheritance *_parent
            ) : name(_name),description(_desc),parentname(_parentname),
                maker(_maker),parent(_parent),
                virt(_maker==0),registered(false) {}


            const STRING& get_name() const {
                return name;
            }


            const STRING& get_description() const {
                return description;
            }


            bool is_virtual() const {
                return virt;
            }


            bool is_instantiable() const {
                return !virt;
            }




            Model* get_model() const {
                return get_model(true);
            }


            Model* get_named_model(
                const STRING& _name
            ) const {
                return get_named_model(_name,false);
            }


            const Inheritance* get_parent() const {
                return parent;
            }



            const InheritanceList& get_children() const {
                return children;
            }



            InheritanceList get_progeny() const {
                InheritanceList result;
                get_model_inheritance_list(result,1);
                return result;
            }


			const Inheritance*
			get_named_inheritance(
				const STRING& model,
				bool nothrow = false
			) const;


            ModelParameterList get_parameters() const {
                ModelParameterList result;
                get_parameters(result,true);
                return result;
            }


            const ModelParameterBase* get_modelparameter(
                const STRING& param
            ) const;


            Model* make() const;


            Model* clone(
                const Model& m
            ) const;



            void Register_Model(Inheritance& topModel = Model::inheritance) {
                if (!registered) {
					if (parent) {


						Inheritance* in = const_cast<Inheritance*>(topModel.get_named_inheritance(parentname));
						if (in) in->add_child(this);
					}
                    registered=true;
                }
            }


        private:
            friend class ModelParameterBase;


            void add_child(
                Inheritance* ib
            ) {
                children.push_back(ib);
            }


            void add_parameter(
                ModelParameterBase* mp
            ) {
                parameters.push_back(mp);
            }


            STRING name;


            STRING description;


			STRING parentname;


            InheritanceList children;


            Inheritance* parent;


            ModelParameterList parameters;


            bool virt;


            bool registered;



            static bool show_model_mode;



            std::unique_ptr<Model_Maker_Base> maker;








            void get_model_inheritance_list(InheritanceList& result,int option=1) const;
            Model* get_model(bool top) const;
            Model* get_named_model(const STRING& name,bool nothrow) const;
			void get_parameters(ModelParameterList& result, bool top) const;
			friend void display_descriptions(const Inheritance* inherit,
                                             const InheritanceList& ilist);
    };







    template <class _Model>
    class Model_Ptr
    {
        public:
            typedef _Model Type;
            typedef _Model* Ptr;


            Model_Ptr() : model(0) {}


            Model_Ptr(const Model_Ptr<_Model>& model_ptr) {
                if (model_ptr.model) {
                    const Inheritance& i = model_ptr.model->get_inheritance();
                    model = static_cast<_Model*>(i.clone(*(model_ptr.model)));
                }
                else  model = 0;
            }


            Model_Ptr& operator=(const Model_Ptr<_Model>& model_ptr) {
                if (model) delete model;
                if (model_ptr.model) {
                    const Inheritance& i = model_ptr.model->get_inheritance();
                    model = static_cast<_Model*>(i.clone(*(model_ptr.model)));
                }
                else model=0;
                return *this;
            }



            Model_Ptr(Model* m) : model(0) {
                *this = m;
            }



            Model_Ptr& operator=(Model* m) {
                if (model) delete model;
                if (m) {

                    model = static_cast<_Model*>(m);
                }
                else model = 0;
                return *this;
            }





            Model_Ptr(const Get_Model_Ptr& gm) : model(0) {
                *this = gm;
            }





            Model_Ptr& operator=(const Get_Model_Ptr&         ) {
                GetPtr();
                return *this;
            }




            Model_Ptr& operator=(const STRING& name) {
                if (model) delete model;
                model = _Get_Named_Model<_Model>(name);
                return *this;
            }




            Model_Ptr(const STRING& name) : model(0) {
                *this = name;
            }



            ~Model_Ptr() {
                if (model) delete model;
            }


            _Model* operator->() const {
                if (!model) throw(SCATMECH_exception(
                                          "Attempt to access element by null pointer"));
                return get();
            }


            _Model& operator*() const {
                return *(get());
            }


            _Model* get() const {
                return model;
            }


            Model_Ptr& GetPtr() {
                if (model) delete model;
                model =  _Get_Model<_Model>();
                return *this;
            }


            static const Inheritance* GetInheritance() {
                return &(_Model::inheritance);
            }

        private:

            _Model* model;
    };






    class Get_Model_Ptr
    {};






    class Model_Maker_Base {
        public:



            virtual Model* make() const = 0;




            virtual Model* clone(const Model& model) const = 0;

            virtual ~Model_Maker_Base() {}
    };




    template <class _Model>
    class Model_Maker : public Model_Maker_Base {
        public:
            Model* make() const {
                Model* newModel = static_cast<Model*>(new _Model);
                newModel->init();
                return newModel;
            }
            Model* clone(const Model& model) const {
                return new _Model(static_cast<const _Model&>(model));
            }

            virtual ~Model_Maker<_Model>() {}
    };




    template <>
    class Model_Maker<void> : public Model_Maker_Base {
        public:
            Model* make() const {
                return static_cast<Model*>(0);
            }
            Model* clone(const Model&            ) const {
                return static_cast<Model*>(0);
            }
            virtual ~Model_Maker<void>() {}
    };

    inline Model* Inheritance::make() const {
        return maker->make();
    }
    inline Model* Inheritance::clone(const Model& m) const {
        return maker->clone(m);
    }






    void Register(
        const Model* x
    );





    template <class _Model>
    _Model* _Get_Model()
    {
        Register((_Model*)0);
        return static_cast<_Model*>(_Model::inheritance.get_model());
    }





    template <class _Model>
    _Model* _Get_Named_Model(
        const STRING& name
    ) {
        Register((_Model*)0);
        return static_cast<_Model*>(_Model::inheritance.get_named_model(name));
    }


#define Register_Model(_Model) _Model::inheritance.Register_Model()





    inline void
    AssertNullParameter(const STRING& parameter)
    {
        if (parameter!="") throw SCATMECH_exception("No subparameter allowed: " + parameter);
    }








    template <class TYPE>
    void ModelParameterSet(
        TYPE& variable,
        const STRING& subparameter,
        const STRING& value
    )
    {
        AssertNullParameter(subparameter);
        variable = from_string<TYPE>(value);
    }








    template <class TYPE>
    STRING ModelParameterGet(
        TYPE& variable,
        const STRING& subparameter
    )
    {
        AssertNullParameter(subparameter);
        return to_string(variable);
    }







    template <class TYPE>
    void ModelParameterAskUser(
        TYPE& variable,
        const STRING& prompt
    )
    {
        variable = AskUser(prompt,variable);
    }




    template <class MODEL,class TYPE>
    class ModelParameter : public ModelParameterBase {
        public:
            ModelParameter(const STRING& _name,
                           const STRING& _description,
                           const STRING& _type,
                           const STRING& _defaultvalue,
                           int _recalclevel,
                           Inheritance& inherit,
                           TYPE MODEL::*_parameter)
                : ModelParameterBase(_name,_description,_type,_defaultvalue,_recalclevel,inherit),
                  parameter(_parameter) {}

            virtual void set_parameter(Model* model,
                                       const STRING& subparameter,
                                       const STRING& value) const
            {
                TYPE &_parameter = ((MODEL*)model)->*parameter;
                ModelParameterSet<TYPE>(_parameter,subparameter,value);
                model->set_recalc(recalclevel);
            }

            virtual void AskUser(Model* model) const {
                TYPE &_parameter = ((MODEL*)model)->*parameter;
                ModelParameterAskUser<TYPE>(_parameter,description);
                model->set_recalc(recalclevel);
            }

            virtual void* get_ptr(const Model* model) const {
                return (void*)&(((MODEL*)model)->*parameter);
            }

            virtual STRING get_parameter(const Model* model,
                                         const STRING& subparameter) const
            {
                TYPE &_parameter = ((MODEL*)model)->*parameter;
                return ModelParameterGet<TYPE>(_parameter,subparameter);
            }

        private:

            TYPE MODEL::*parameter;
    };



    template <class MODEL,class TYPE>
    class ModelPtrParameter : public ModelParameterBase {
        public:
            ModelPtrParameter(const STRING& _name,
                              const STRING& _description,
                              const STRING& _type,
                              const STRING& _defaultvalue,
                              int _recalclevel,
                              Inheritance& inherit,
                              TYPE MODEL::*_parameter)
                : ModelParameterBase(_name,_description,_type,_defaultvalue,_recalclevel,inherit),
                  parameter(_parameter) {}

            virtual void set_parameter(Model* model,
                                       const STRING& subparameter,
                                       const STRING& value) const
            {
                TYPE &_parameter = ((MODEL*)model)->*parameter;

                if (subparameter=="") {


                    Register(TYPE().get());
                    _parameter = TYPE::GetInheritance()->get_named_model(value);
                } else {

                    _parameter->set_parameter(subparameter,value);
                }
                model->set_recalc(recalclevel);
            }

            virtual void AskUser(Model* model) const {
                TYPE &_parameter = ((MODEL*)model)->*parameter;
                SCATMECH_output << description;
                _parameter = Get_Model_Ptr();
				_parameter->Model::AskUser();
                model->set_recalc(recalclevel);
            }

            virtual const Inheritance* get_inheritance(const Model* model) const {
                if (model!=NULL) {
                    return &((((MODEL*)model)->*parameter).get()->get_inheritance());
                } else {
                    return TYPE::GetInheritance();
                }
            }

            virtual void* get_ptr(const Model* model) const {
                return (void*)((((MODEL*)model)->*parameter).get());
            }

            virtual STRING get_parameter(const Model* model,
                                         const STRING& subparameter) const
            {
                if (subparameter=="") {
                    return ((MODEL*)model->*parameter).get()->get_inheritance().get_name();
                } else {
                    return ((MODEL*)model->*parameter).get()->get_parameter(subparameter);
                }
            }

        private:


            TYPE MODEL::*parameter;
	};








#define DECLARE_MODEL() \
		public: static SCATMECH::Inheritance inheritance; \
		public: SCATMECH::Model* clone() const {return inheritance.clone(*this);} \
		public: virtual const SCATMECH::Inheritance& get_inheritance() const \
						{return inheritance;}











#define DECLARE_PARAMETER(TYPE,PARAMETER) \
		public:    const TYPE& get_##PARAMETER() const {return PARAMETER;} \
		public:    void set_##PARAMETER(const TYPE& __x) {PARAMETER = __x; set_recalc(MP_##PARAMETER->recalclevel);} \
		protected: TYPE PARAMETER; \
		private: static SCATMECH::ModelParameterBasePtr MP_##PARAMETER;












#define DEFINE_MODEL(child,parent,desc) \
	SCATMECH::Inheritance child::inheritance(#child, \
		desc, #parent, \
		new SCATMECH::Model_Maker<child>(), \
		&parent::inheritance);












#define DEFINE_VIRTUAL_MODEL(child,parent,desc) \
	SCATMECH::Inheritance child::inheritance(#child, \
		desc, #parent, \
		0, \
		&parent::inheritance);















#define DEFINE_PARAMETER(CLASS,TYPE,PARAMETER,DESCRIPTION,DEFAULTVALUE,RECALCLEVEL) \
		ModelParameterBasePtr \
		CLASS::MP_##PARAMETER( \
							  new ModelParameter< CLASS , TYPE >(#PARAMETER, \
																 DESCRIPTION, \
																 #TYPE, \
																 DEFAULTVALUE, \
																 RECALCLEVEL, \
																 CLASS::inheritance, \
																 &CLASS::PARAMETER \
																) \
							  );















#define DEFINE_PTRPARAMETER(CLASS,TYPE,PARAMETER,DESCRIPTION,DEFAULTVALUE,RECALCLEVEL) \
		ModelParameterBasePtr \
		CLASS::MP_##PARAMETER( \
							  new ModelPtrParameter< CLASS , TYPE >(#PARAMETER, \
																	DESCRIPTION, \
																	#TYPE, \
																	DEFAULTVALUE, \
																	RECALCLEVEL, \
																	CLASS::inheritance, \
																	&CLASS::PARAMETER \
																   ) \
							 );

}

#endif
