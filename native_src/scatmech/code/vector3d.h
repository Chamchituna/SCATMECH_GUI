











#ifndef SCATMECH_VECTOR3D_H
#define SCATMECH_VECTOR3D_H

#include "scatmech.h"
#include <ios>

namespace SCATMECH {




    template <class TYPE>
    class Vector3D {
        public:



            TYPE x,y,z;


            Vector3D<TYPE>() {}


            Vector3D<TYPE>(const TYPE& _x,const TYPE& _y,const TYPE& _z)
            {
                x=_x;
                y=_y;
                z=_z;
            }


            template <class TYPE2>
            Vector3D<TYPE>(const Vector3D<TYPE2>& a)
            {
                x=a.x;
                y=a.y;
                z=a.z;
            }


            Vector3D<TYPE>& operator=(const Vector3D<TYPE>& a)
            {
                x=a.x;
                y=a.y;
                z=a.z;
                return *this;
            }


            Vector3D<TYPE> operator+(const Vector3D<TYPE>& a) const
            {
                return Vector3D<TYPE>(x+a.x,y+a.y,z+a.z);
            }


            const Vector3D<TYPE>& operator+=(const Vector3D<TYPE>& a)
            {
                return *this = *this+a;
            }


            Vector3D<TYPE> operator-(const Vector3D<TYPE>& a) const
            {
                return Vector3D<TYPE>(x-a.x,y-a.y,z-a.z);
            }


            const Vector3D<TYPE>& operator-=(const Vector3D<TYPE>& a)
            {
                return *this = *this-a;
            }


            Vector3D<TYPE> operator-() const
            {
                return Vector3D<TYPE>(-x,-y,-z);
            }


            TYPE operator*(const Vector3D<TYPE>& a) const
            {
                return x*a.x+y*a.y+z*a.z;
            }


            Vector3D<TYPE> operator*(const TYPE& b) const
            {
                return Vector3D<TYPE>(x*b,y*b,z*b);
            }


            friend Vector3D<TYPE> operator*(const TYPE& b,
                                            const Vector3D<TYPE>& a)
            {
                return Vector3D<TYPE>(a.x*b,a.y*b,a.z*b);
            }


            Vector3D<TYPE> operator/(const TYPE& a) const
            {
                return Vector3D<TYPE>(x/a,y/a,z/a);
            }


            const Vector3D<TYPE>& operator/=(const TYPE& a)
            {
                return *this = *this/a;
            }


            friend Vector3D<TYPE> cross(const Vector3D<TYPE>& a,
                                        const Vector3D<TYPE>& b)
            {   return Vector3D<TYPE>(a.y*b.z-a.z*b.y,
                                      a.z*b.x-a.x*b.z,
                                      a.x*b.y-a.y*b.x);
            }


            const Vector3D<TYPE>& operator*=(const TYPE& a)
            {
                return *this = *this*a;
            }


            bool operator==(const Vector3D<TYPE>& b) const
            {
                return (x==b.x && y==b.y && z==b.z);
            }


            bool operator!=(const Vector3D<TYPE>& b) const
            {
                return (x!=b.x || y!=b.y || z!=b.z);
            }


            friend std::ostream& operator<<(std::ostream& os,const Vector3D<TYPE>& v)
            {
                return os << '(' << v.x << ',' << v.y << ',' << v.z << ')';
            }


            friend std::istream& operator>>(std::istream& is,Vector3D<TYPE>& v)
            {   TYPE vx,vy,vz;
                char c;
                is >> c;
                if (c=='(') {
                    is >> vx;
                    is >> c;
                    if (c==',') {
                        is >> vy;
                        is >> c;
                        if (c==',') {
                            is >> vz;
                            is >> c;
                            if (c==')') {
                                v.x=vx, v.y=vy, v.z=vz;
                                return is;
                            }
                        }
                    }
                }

                is.setstate(std::ios::failbit);
                return is;
            }
            const TYPE& operator[](int i) const {
                return *(&x+i);
            }
            TYPE& operator[](int i) {
                return *(&x+i);
            }

    };


    inline
    double
    Norm(const Vector3D<double>& a)
    {
        return sqrt(sqr(a.x)+sqr(a.y)+sqr(a.z));
    }


    inline
    Vector3D<double>
    unit(const Vector3D<double>& a)
    {
        return a/Norm(a);

    }


    inline
    Vector3D<double>
    polar(double r,double theta,double phi)
    {
        return Vector3D<double>(r*sin(theta)*cos(phi),
                                r*sin(theta)*sin(phi),
                                r*cos(theta));
    }


    inline
    Vector3D<COMPLEX>
    Conj(const Vector3D<COMPLEX>& a)
    {
        using std::conj;
        return Vector3D<COMPLEX>(conj(a.x),conj(a.y),conj(a.z));
    }


    inline
    double
    Norm(const Vector3D<COMPLEX >& a)
    {
        using std::norm;
        return sqrt(norm(a.x)+norm(a.y)+norm(a.z));
    }


    inline
    Vector3D<COMPLEX >
    unit(const Vector3D<COMPLEX >& a)
    {
        return a/Norm(a);
    }


    Vector3D<double> perpto(const Vector3D<double>& a,
                            const Vector3D<double>& b);
    Vector3D<COMPLEX > perpto(const Vector3D<COMPLEX >& a,
                              const Vector3D<COMPLEX >& b);


    typedef Vector3D<double> Vector;

    typedef Vector3D<COMPLEX > CVector;


}

#endif
