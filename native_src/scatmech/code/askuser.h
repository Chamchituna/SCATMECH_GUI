











#ifndef SCATMECH_ASKUSER_H
#define SCATMECH_ASKUSER_H

#include <fstream>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <stack>
#include <list>
#include <cctype>
#include <cstddef>
#include <algorithm>

namespace SCATMECH {

    extern std::string SCATMECH_DIRECTORY;
    void check_SCATMECH_DIRECTORY();
    std::string find_file(const std::string& name);







    class streambuf_with_comments : public std::streambuf {
        private:
            std::streambuf* istr;
            char comment_delimiter;

        protected:
            void clear_comment();
            int underflow();
            int uflow();
            int pbackfail(char c);

            std::streampos seekpos(std::streampos sp, std::ios_base::openmode which) {
                return istr->pubseekpos(sp,which);
            }
            std::streampos seekoff(std::streamoff off, std::ios_base::seekdir way,std::ios_base::openmode which) {
                return istr->pubseekoff(off,way,which);
            }
            int sync() {
                return istr->pubsync();
            }

        public:
            streambuf_with_comments(std::streambuf *s);
            streambuf_with_comments& operator=(std::streambuf *s);
            virtual ~streambuf_with_comments();
            void set_comment_delimiter(char c) {
                comment_delimiter=c;
            }
    };






    class streambuf_with_echo : public std::streambuf {
        private:
            std::streambuf* strbuf;
            std::ofstream echostream;
            bool echo;

        protected:
            int underflow();
            int uflow();
            int pbackfail(char c);

        public:
            streambuf_with_echo(std::streambuf *s,const std::string& name="");
            streambuf_with_echo& operator=(std::streambuf *s);
            virtual ~streambuf_with_echo();
            void set_echo(const std::string& name);
            void flush_echo() {
                echostream.flush();
            }
    };






    class istream_with_comments : public std::istream {
        public:
            typedef std::istream istream;
            istream_with_comments(std::streambuf* _buf) : buf(_buf), istream(&buf) {}

            istream_with_comments& operator=(std::istream& _str) {
                buf = _str.rdbuf();
                return *this;
            }


            std::string getstrline();


            std::string getquoted();

        protected:
            streambuf_with_comments buf;
    };






    class ifstream_with_comments : public istream_with_comments {
        public:
            ifstream_with_comments(const char* name) : fbuf(), istream_with_comments(&fbuf) {
                using namespace std;
                std::string newname = find_file(name);
                open(newname.c_str());
            }
            ifstream_with_comments() : istream_with_comments(&fbuf) {}
            void open(const char* name) {
                using namespace std;
                fbuf.open(name,ios::in);
                if (!fbuf.is_open()) setstate(failbit);
            }
            void close() {
                fbuf.close();
            }
        protected:
            std::filebuf fbuf;
    };


    extern istream_with_comments SCATMECH_input;


    extern std::ostream SCATMECH_output;






    template <class TYPE>
    TYPE AskUser(
        const std::string& prompt,
        const TYPE& deflt
    )
    {
        while (1) {
            try {
                SCATMECH_output << prompt << " <" << deflt << "> :";

                std::istringstream line(SCATMECH_input.getstrline());

                TYPE result;
                line >> result;

                if (line.fail()) result = deflt;

                return result;
            }
            catch (std::exception& e) {
                SCATMECH_output << e.what() << std::endl;
                if (SCATMECH_input.eof()) throw e;
            }
        }
    }


    inline std::string AskUser(const std::string& query, const char* deflt)
    {
        return AskUser(query,std::string(deflt));
    }


    void set_istream(std::istream &is);


    void set_ostream(std::ostream &os);






    template <class TYPE>
    std::string format(const std::string& fmt, const TYPE& value) {
        char buffer[256];
        #ifdef _MSC_VER
        sprintf_s(buffer,256,fmt.c_str(),value);
        #else
        snprintf(buffer,256,fmt.c_str(),value);
        #endif
        return std::string(buffer);
    }


    std::string getstrline(std::istream& str);




	class CommandLineParser {
	private:
		typedef std::list<std::string> TokenList;
	public:
		CommandLineParser(const int &argc, const char **argv) {
			for (int i = 1; i < argc; ++i)
				tokens.push_back(argv[i]);
		}
		std::string getOption(const std::string &option) const {
			TokenList::const_iterator itr;
			itr = std::find(tokens.begin(), tokens.end(), option);
			if (itr != tokens.end() && ++itr != tokens.end()) {
				return *itr;
			}
			else {
				return std::string("");
			}
		}
		bool optionExists(const std::string &option) const {
			return (std::find(tokens.begin(), tokens.end(), option)!= tokens.end());
		}
	private:
		TokenList tokens;
	};



}


#endif
