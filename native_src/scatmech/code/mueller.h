











#ifndef SCATMECH_MUELLER_HPP
#define SCATMECH_MUELLER_HPP

#include "scatmech.h"
#include "optconst.h"
#include "vector3d.h"












namespace SCATMECH {



    class JonesMatrix;
    class JonesVector;
    class MuellerMatrix;
    class StokesVector;







    class JonesVector {
        public:




            JonesVector() {}
            JonesVector(const JonesVector& x)
            {
                j[0]=x.j[0];
                j[1]=x.j[1];
            }



            JonesVector(const COMPLEX& s,const COMPLEX& p)
            {
                j[0]=s;
                j[1]=p;
            }





            explicit JonesVector(const StokesVector& x);




            JonesVector& operator=(const JonesVector& x)
            {
                j[0]=x.j[0];
                j[1]=x.j[1];
                return *this;
            }




            JonesVector operator*(const COMPLEX& x) const
            {
                return JonesVector(j[0]*x,j[1]*x);
            }
            friend JonesVector operator*(const COMPLEX& x,const JonesVector& y)
            {
                return y*x;
            }
            JonesVector& operator*=(const COMPLEX& x)
            {
                return *this = x*(*this);
            }



            JonesVector operator/(const COMPLEX& x) const
            {
                return (1./x)*(*this);
            }
            JonesVector& operator/=(const COMPLEX& x)
            {
                return *this = (1./x)*(*this);
            }





            COMPLEX operator*(const JonesVector& v) const
            {
                return j[0]*v[0]+j[1]*v[1];
            }



            friend JonesVector conj(const JonesVector& v)
            {
                return JonesVector(conj(v[0]),conj(v[1]));
            }




            JonesVector operator+(const JonesVector& a) const
            {
                return JonesVector(j[0]+a.j[0],j[1]+a.j[1]);
            }
            JonesVector operator-(const JonesVector& a) const
            {
                return JonesVector(j[0]-a.j[0],j[1]-a.j[1]);
            }
            JonesVector operator+=(const JonesVector& a)
            {
                return *this = *this + a;
            }
            JonesVector operator-=(const JonesVector& a)
            {
                return *this = *this - a;
            }
            JonesVector operator-() const {
                return JonesVector(-j[0],-j[1]);
            }
            const JonesVector& operator+() const {
                return *this;
            }




            COMPLEX& operator[](int i) {
                return j[i];
            }
            const COMPLEX& operator[](int i) const {
                return j[i];
            }




            double intensity() const {
                return std::norm(j[0])+std::norm(j[1]);
            }




            double psi() const
            {
                return atan(std::abs(j[1])/std::abs(j[0]));
            }




            double delta() const
            {
                return std::arg(j[1])-std::arg(j[0]);
            }




            double eta() const;




            double DOLP() const;




            double DOP() const
            {
                return 1.;
            }




            double DOCP() const;


            double e() const;




            double epsilon() const
            {
                double t = e();
                return sqrt(1.-t*t);
            }




            JonesVector rotate(const double angle) const;




            COMPLEX& S() {
                return j[0];
            }
            const COMPLEX& S() const {
                return j[0];
            }




            COMPLEX& P() {
                return j[1];
            }
            const COMPLEX& P() const {
                return j[1];
            }

            friend class JonesMatrix;
            friend class MuellerMatrix;
            friend class StokesVector;

        private:



            COMPLEX j[2];
    };







    class JonesMatrix {
        public:



            JonesMatrix() {}




            JonesMatrix(const JonesMatrix& x)
            {
                j[0]=x.j[0];
                j[1]=x.j[1];
                j[2]=x.j[2];
                j[3]=x.j[3];
            }




            JonesMatrix(const COMPLEX& pp,const COMPLEX& ss,
                        const COMPLEX& ps,const COMPLEX& sp)
            {
                j[0]=pp;
                j[1]=ss;
                j[2]=ps,
                j[3]=sp;
            }








            explicit JonesMatrix(const MuellerMatrix& x);




            JonesMatrix& operator=(const JonesMatrix& x)
            {
                j[0]=x.j[0];
                j[1]=x.j[1];
                j[2]=x.j[2];
                j[3]=x.j[3];
                return *this;
            }




            JonesMatrix operator*(const JonesMatrix& matrix) const
            {
                return JonesMatrix(j[3]*matrix.j[2]+j[0]*matrix.j[0],
                                   j[1]*matrix.j[1]+j[2]*matrix.j[3],
                                   j[1]*matrix.j[2]+j[2]*matrix.j[0],
                                   j[3]*matrix.j[1]+j[0]*matrix.j[3]);
            }
            JonesMatrix& operator*=(const JonesMatrix& a)
            {
                return *this = *this * a;
            }




            JonesMatrix operator*(const COMPLEX& x) const
            {
                return JonesMatrix(j[0]*x,j[1]*x,j[2]*x,j[3]*x);
            }
            friend JonesMatrix operator*(const COMPLEX& x,const JonesMatrix& y)
            {
                return y*x;
            }
            JonesMatrix& operator*=(const COMPLEX& a) {
                return *this = *this * a;
            }




            JonesMatrix operator/(const COMPLEX& x) const
            {
                return JonesMatrix(j[0]/x,j[1]/x,j[2]/x,j[3]/x);
            }
            JonesMatrix& operator/=(const COMPLEX& a) {
                return *this = *this / a;
            }




            JonesVector operator*(const JonesVector& a) const
            {
                return JonesVector(j[1]*a.j[0]+j[2]*a.j[1],j[3]*a.j[0]+j[0]*a.j[1]);
            }




            JonesMatrix operator+(const JonesMatrix& a) const
            {
                return JonesMatrix(j[0]+a.j[0],j[1]+a.j[1],j[2]+a.j[2],j[3]+a.j[3]);
            }
            JonesMatrix operator-(const JonesMatrix& a) const
            {
                return JonesMatrix(j[0]-a.j[0],j[1]-a.j[1],j[2]-a.j[2],j[3]-a.j[3]);
            }
            JonesMatrix& operator+=(const JonesMatrix& a)
            {
                return *this = *this + a;
            }
            JonesMatrix& operator-=(const JonesMatrix& a)
            {
                return *this = *this - a;
            }
            JonesMatrix operator-() const {
                return JonesMatrix(-j[0],-j[1],-j[2],-j[3]);
            }
            const JonesMatrix& operator+() const {
                return *this;
            }




            COMPLEX& operator[](int i) {
                return j[i];
            };
            const COMPLEX& operator[](int i) const {
                return j[i];
            };




            JonesMatrix rotate(double angle) const;




            JonesMatrix transpose() const {
                return JonesMatrix(j[0],j[1],j[3],j[2]);
            }
            friend JonesMatrix transpose(const JonesMatrix& j) {
                return JonesMatrix(j[0],j[1],j[3],j[2]);
            }




            JonesMatrix hermitian() const {
                return JonesMatrix(std::conj(j[0]), std::conj(j[1]), std::conj(j[3]), std::conj(j[2]));
            }
            friend JonesMatrix hermitian(const JonesMatrix& j) {
                return JonesMatrix(std::conj(j[0]), std::conj(j[1]), std::conj(j[3]), std::conj(j[2]));
            }




            JonesMatrix conj() const {
                return JonesMatrix(std::conj(j[0]), std::conj(j[1]), std::conj(j[2]), std::conj(j[3]));
            }
            friend JonesMatrix conj(const JonesMatrix& j) {
                return JonesMatrix(std::conj(j[0]), std::conj(j[1]), std::conj(j[2]), std::conj(j[3]));
            }




            COMPLEX& PP() {
                return j[0];
            }
            const COMPLEX& PP() const {
                return j[0];
            }




            COMPLEX& SS() {
                return j[1];
            }
            const COMPLEX& SS() const {
                return j[1];
            }




            COMPLEX& PS() {
                return j[2];
            }
            const COMPLEX& PS() const {
                return j[2];
            }




            COMPLEX& SP() {
                return j[3];
            }
            const COMPLEX& SP() const {
                return j[3];
            }


            friend class JonesVector;
            friend class MuellerMatrix;
            friend class StokesVector;
        private:




            COMPLEX j[4];
    };






    class MuellerMatrix {

        public:




            MuellerMatrix() {}




            MuellerMatrix(const JonesMatrix& j);




            MuellerMatrix(const MuellerMatrix& x) {
                *this=x;
            }






            MuellerMatrix(const StokesVector& s1,const StokesVector& s2,const StokesVector& s3, const StokesVector& s4);




            MuellerMatrix& operator=(const MuellerMatrix& x);




            double* operator[](int i) {
                return m[i];
            }
            const double* operator[](int i) const {
                return m[i];
            }




            double Tmax() const
            {
                return m[0][0]+sqrt(m[0][1]*m[0][1]+m[0][2]*m[0][2]+m[0][3]*m[0][3]);
            }



            double Tmin() const
            {
                return m[0][0]-sqrt(m[0][1]*m[0][1]+m[0][2]*m[0][2]+m[0][3]*m[0][3]);
            }



            double diattenuation() const
            {
                return sqrt(m[0][1]*m[0][1]+m[0][2]*m[0][2]+m[0][3]*m[0][3])/m[0][0];
            }



            double linear_diattenuation() const
            {
                return sqrt(m[0][1]*m[0][1]+m[0][2]*m[0][2])/m[0][0];
            }



            double polarization_dependent_loss() const
            {
                return 10.*std::log(Tmax()/Tmin())/std::log(10.);
            }



            double polarizance() const
            {
				return sqrt(m[1][0] * m[1][0] + m[2][0] * m[2][0] + m[3][0] * m[3][0]) / m[0][0];
            }



            double depolarization_index() const
            {
                MuellerMatrix mmt = (*this)*(this->transpose());
                return sqrt((mmt[0][0]+mmt[1][1]+mmt[2][2]+mmt[3][3]-sqr(m[0][0]))/(3.*sqr(m[0][0])));
            }



            double extinction_ratio() const {
                return Tmax()/Tmin();
            }




            MuellerMatrix rotate(const double angle) const;




            MuellerMatrix parity() const;




            MuellerMatrix transpose() const;




            MuellerMatrix inverse() const;




            bool valid() const;








            bool physically_valid() const;




            MuellerMatrix log() const;




            MuellerMatrix exp() const;




            void Cloude_Decomposition(MuellerMatrix& M1,MuellerMatrix& M2,MuellerMatrix& M3,MuellerMatrix& M4) const;




            void Lu_Chipman_Decomposition(MuellerMatrix& depolarizer,MuellerMatrix& retarder, MuellerMatrix& diattenuator) const;




            MuellerMatrix Closest_NonDepolarizing(int rank=1) const;




            MuellerMatrix normalized() const {
                return *this/m[0][0];
            }




			double entropy() const;




            MuellerMatrix operator*(const MuellerMatrix& matrix) const;
            MuellerMatrix& operator*=(const MuellerMatrix& matrix)
            {
                return *this = *this * matrix;
            }




            StokesVector operator*(const StokesVector &s) const;




            MuellerMatrix operator*(double d) const;
            friend MuellerMatrix operator*(double d,const MuellerMatrix &v)
            {
                return v*d;
            }
            MuellerMatrix& operator*=(const double d) {
                return *this = *this * d;
            }




            MuellerMatrix operator/(double d) const;
            MuellerMatrix& operator/=(double d) {
                return *this = *this * (1./d);
            }




            MuellerMatrix operator+(const MuellerMatrix& a) const;
            MuellerMatrix& operator+=(const MuellerMatrix& matrix)
            {
                return *this = *this + matrix;
            }




            MuellerMatrix operator-(const MuellerMatrix& a) const;
            MuellerMatrix& operator-=(const MuellerMatrix& matrix)
            {
                return *this = *this - matrix;
            }

            MuellerMatrix operator-() const;
            const MuellerMatrix& operator+() const {
                return *this;
            }

            friend class JonesVector;
            friend class JonesMatrix;
            friend class StokesVector;

        private:




            double m[4][4];
    };






    class StokesVector {

        public:




            StokesVector() {}




            StokesVector(double I,double Q,double U,double V)
            {
                s[0]=I;
                s[1]=Q;
                s[2]=U;
                s[3]=V;
            }




            StokesVector(const StokesVector& x)
            {
                s[0]=x.s[0];
                s[1]=x.s[1];
                s[2]=x.s[2];
                s[3]=x.s[3];
            }




            StokesVector(const JonesVector& j)
            {
                s[0] = std::norm(j.j[0])+std::norm(j.j[1]);
                s[1] = std::norm(j.j[0])-std::norm(j.j[1]);
                s[2] = 2*std::real(std::conj(j.j[0])*j.j[1]);
                s[3] = 2*std::imag(std::conj(j.j[0])*j.j[1]);
            }




            StokesVector& operator=(const StokesVector& x)
            {
                s[0]=x.s[0];
                s[1]=x.s[1];
                s[2]=x.s[2];
                s[3]=x.s[3];
                return *this;
            }




            double& I() {
                return s[0];
            }
            double I() const {
                return s[0];
            }




            double& Q() {
                return s[1];
            }
            double Q() const {
                return s[1];
            }




            double& U() {
                return s[2];
            }
            double U() const {
                return s[2];
            }




            double& V() {
                return s[3];
            }
            double V() const {
                return s[3];
            }




            double q() const {
                return Q()/I();
            }




            double u() const {
                return U()/I();
            }




            double v() const {
                return V()/I();
            }




            double& operator[](int i) {
                return s[i];
            }
            double operator[](int i) const {
                return s[i];
            }




            StokesVector operator*(const MuellerMatrix& matrix) const;




            double operator*(const StokesVector& a) const
            {
                return a.s[0]*s[0] + a.s[1]*s[1] + a.s[2]*s[2] + a.s[3]*s[3];
            }




            StokesVector operator+(const StokesVector& a) const
            {
                return StokesVector(s[0]+a.s[0],s[1]+a.s[1],s[2]+a.s[2],s[3]+a.s[3]);
            }
            StokesVector operator+=(const StokesVector& a)
            {
                return *this = *this+a;
            }




            StokesVector operator-(const StokesVector& a) const
            {
                return StokesVector(s[0]-a.s[0],s[1]-a.s[1],s[2]-a.s[2],s[3]-a.s[3]);
            }
            StokesVector operator-=(const StokesVector& a)
            {
                return *this = *this-a;
            }




            StokesVector operator-() const
            {
                return StokesVector(-s[0],-s[1],-s[2],-s[3]);
            }
            const StokesVector& operator+() const {
                return *this;
            }




            StokesVector operator*(double d) const
            {
                return StokesVector(s[0]*d,s[1]*d,s[2]*d,s[3]*d);
            }
            friend StokesVector operator*(double d,const StokesVector& s)
            {
                return StokesVector(s[0]*d,s[1]*d,s[2]*d,s[3]*d);
            }
            StokesVector& operator*=(double d) {
                return *this = *this *  d;
            }




            StokesVector operator/(double d) const
            {
                return StokesVector(s[0]/d,s[1]/d,s[2]/d,s[3]/d);
            }
            StokesVector& operator/=(double d) {
                return *this = *this /  d;
            }




            StokesVector rotate(double angle) const;




            double eta() const
            {
                return atan2(s[2],s[1])/2.;
            }




            double intensity() const
            {
                return s[0];
            }




            double DOLP() const
            {
                return sqrt(sqr(s[2])+sqr(s[1]))/s[0];
            }




            double DOP() const
            {
                return sqrt(sqr(s[2])+sqr(s[1])+sqr(s[3]))/s[0];
            }




            double DOCP() const
            {
                return s[3]/s[0];
            }




            double e() const
            {
                return s[3]/(s[0]+sqrt(sqr(s[1])+sqr(s[2])));
            }




            double delta() const;




            double psi() const;




            double eccentricity() const
            {
                return sqrt(1.-sqr(e()));
            }




            bool valid() const
            {
                return s[0]>=sqrt(sqr(s[1])+sqr(s[2])+sqr(s[3]));
            }




            StokesVector pol_part() const;




            StokesVector unpol_part() const;

            friend class JonesVector;
            friend class MuellerMatrix;
            friend class JonesMatrix;

        private:




            double s[4];
    };




    std::ostream& operator<<(std::ostream& os,const JonesVector& j);
    std::istream& operator>>(std::istream& is,JonesVector& j);
    std::ostream& operator<<(std::ostream& os,const StokesVector& j);
    std::istream& operator>>(std::istream& is, StokesVector& j);
    std::ostream& operator<<(std::ostream& os,const MuellerMatrix& mm);
    std::istream& operator>>(std::istream& is, MuellerMatrix& mm);
    std::ostream& operator<<(std::ostream& os,const JonesMatrix& j);
    std::istream& operator>>(std::istream& is, JonesMatrix& j);




    inline JonesMatrix JonesZero() {
        return JonesMatrix(0.,0.,0.,0.);
    }



    inline JonesMatrix JonesUnit() {
        return JonesMatrix(1.,1.,0.,0.);
    }




    JonesMatrix JonesRotator(double angle);



    JonesMatrix JonesLinearRetarder(double phase, double angle=0);



    JonesMatrix JonesCircularRetarder(double phase);



    JonesMatrix JonesLinearPolarizer(double angle=0,double diattenuation=1);



    JonesMatrix JonesCircularPolarizer(double diattenuation=1);




    JonesMatrix JonesGeneralized(
        const JonesVector& a,
        const JonesVector& b,
        const COMPLEX& ma,
        const COMPLEX& mb
    );




    MuellerMatrix MuellerZero();



    MuellerMatrix MuellerUnit(
        double attenuation=1.
    );



    MuellerMatrix MuellerDepolarizer(
        double attenuation=1.,
        double depolarization=1.
    );



	MuellerMatrix MuellerDiagonal(double m00,double m11, double m22, double m33);



    MuellerMatrix MuellerPartialLinearPolarizer(
        double tmax,
        double tmin,
        double angle
    );




    inline StokesVector StokesZero() {
        return StokesVector(0,0,0,0);
    }




    inline StokesVector StokesVectorUnitUnpol() {
        return StokesVector(1,0,0,0);
    }



    inline StokesVector StokesVectorUnitS() {
        return StokesVector(1,1,0,0);
    }



    inline StokesVector StokesVectorUnitP() {
        return StokesVector(1,-1,0,0);
    }



    inline StokesVector StokesVectorUnitLinear(double eta) {
        return StokesVector(1,cos(2*eta),sin(2*eta),0);
    }



    inline StokesVector StokesVectorUnitRCP() {
        return StokesVector(1,0,0,-1);
    }



    inline StokesVector StokesVectorUnitLCP() {
        return StokesVector(1,0,0,1);
    }



    inline StokesVector StokesVectorUnitGeneral(
        double eta,
        double DOCP=0.,
        double DOP=1.
    )
    {
        double temp=sqrt(1.-sqr(DOCP));
        return StokesVector(1,
                            DOP*temp*cos(2*eta),
                            DOP*temp*sin(2*eta),
                            DOP*DOCP);
    }




    inline JonesVector JonesVectorZero() {
        return JonesVector(0.,0.);
    }



    inline JonesVector JonesVectorUnitS() {
        return JonesVector(1.,0.);
    }



    inline JonesVector JonesVectorUnitP() {
        return JonesVector(0.,1.);
    }



    inline JonesVector JonesVectorUnitLinear(double eta) {
        return JonesVector(cos(eta),sin(eta));
    }



    inline JonesVector JonesVectorUnitRCP() {
        return JonesVector(1.,COMPLEX(0,-1));
    }



    inline JonesVector JonesVectorUnitLCP() {
        return JonesVector(1.,COMPLEX(0,1));
    }



    inline JonesVector JonesVectorUnitGeneral(
        double eta,
        double DOCP
    )
    {
        return JonesVector(StokesVectorUnitGeneral(eta,DOCP));
    }





    void GetBasisVectorsSP(const Vector& k,Vector& s,Vector& p);






    void GetBasisVectorsXY(const Vector& k,Vector& x,Vector& y);






    void GetBasisVectorsParPerp(const Vector& kin,const Vector& kout, Vector& perp, Vector& parin, Vector& parout);



    JonesMatrix GetJonesRotator(const Vector& xo, const Vector& yo, const Vector& xi, const Vector& yi);














    MuellerMatrix ReCrossMueller(const JonesMatrix& j1,const JonesMatrix& j2);





    MuellerMatrix ImCrossMueller(const JonesMatrix& j1,const JonesMatrix& j2);


}

#endif
