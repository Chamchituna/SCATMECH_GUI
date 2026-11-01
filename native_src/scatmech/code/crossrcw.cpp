











#include "crossrcw.h"
#include "matrix3d.h"

using namespace std;

namespace SCATMECH {

#define HERE    message(string("Location: ") + string(__FILE__) + string("  ") + to_string(__LINE__) + "..........");
#define SHOW(file,x)   (file) << #x << " = " << (x) << endl;









    template <class T> void SWAP(T& x,T& y) {
        T t=x;
        x=y;
        y=t;
    }

    static COMPLEX cI(0.,1.);

    void CrossRCW_Model::setup()
    {
        Model::setup();

        if (type<0 || type>3) error("type must be 0, 1, 2, or 3");




        if (RECALC == 0x02) {
            if ((type & 0x02) == (old_type & 0x02)) {
                old_type = type;
                return;
            }
        }

        int i,j,k,l,q;



        if (grating->get_lambda()!=lambda) grating->set_lambda(lambda);
        if (grating->get_order1()!=order1) grating->set_order1(order1);
        if (grating->get_order2()!=order2) grating->set_order2(order2);


        k0 = 2*pi/lambda;
        double mu = 1.;



        CFARRAY E0 = grating->get_EPS0();
        CFARRAY E11 = grating->get_EPS11();
		CFARRAY E12 = grating->get_EPS12();
		CFARRAY E2 = grating->get_EPS2();
        CFARRAY E3 = grating->get_EPS3();

		CFARRAY MU0 = grating->get_MU0();
		CFARRAY MU11 = grating->get_MU11();
		CFARRAY MU12 = grating->get_MU12();
		CFARRAY MU2 = grating->get_MU2();
		CFARRAY MU3 = grating->get_MU3();


        int    levels = grating->get_levels();
        double d1     = grating->get_d1();
        double d2     = grating->get_d2();
        double zeta   = grating->get_zeta()*deg;


        DFARRAY thick(levels,1);
        totalthickness=0;
        for (i=1; i<=levels; ++i) {
            thick(i)=grating->get_thick(i);
            totalthickness+=thick(i);
        }


        double K1 = 2*pi/d1;
        double K2 = 2*pi/d2;

        int M1 = 2*order1+1;
        int M2 = 2*order2+1;
        int MM = M1*M2;
        int MM2 = MM*2;

        E0.array(M1,M2,M1,M2,levels);
        E11.array(M1,M2,M1,M2,levels);
		E12.array(M1, M2, M1, M2, levels);
		E2.array(M1,M2,M1,M2,levels);
        E3.array(M1,M2,M1,M2,levels);
		MU0.array(M1, M2, M1, M2, levels);
		MU11.array(M1, M2, M1, M2, levels);
		MU12.array(M1, M2, M1, M2, levels);
		MU2.array(M1, M2, M1, M2, levels);
		MU3.array(M1, M2, M1, M2, levels);

        double phi = -rotation*deg;
        double cosphi = cos(phi);
        double sinphi = sin(phi);

        double coszeta = cos(zeta);
        double seczeta = 1/coszeta;
        double sinzeta = sin(zeta);
        double tanzeta = tan(zeta);
        double sqrcoszeta = sqr(coszeta);
        double sqrsinzeta = sqr(sinzeta);



        double ep1;
        COMPLEX em1;
        if (type == 0 || type == 1) {
            ep1 = grating->get_medium_i().e1(lambda);
            if (grating->get_medium_i().e2(lambda)!=0) error("Incident medium cannot be absorbing");
            em1 = COMPLEX(grating->get_medium_t().epsilon(lambda));
        } else {
            ep1 = grating->get_medium_t().e1(lambda);
            if (grating->get_medium_t().e2(lambda)!=0) error("Incident medium cannot be absorbing");
            em1 = COMPLEX(grating->get_medium_i().epsilon(lambda));
        }

        eI = ep1;
        eII = em1;
        nI = sqrt(eI);
        nII = sqrt(eII);

        double np1 = sqrt(ep1);
        COMPLEX nm1 = sqrt(em1);

        double kp1 = k0*sqrt(ep1*mu);
        COMPLEX km1 = k0*sqrt(em1*mu);



        double alpha0 = kp1*sin(thetai*deg)*cos(phi);
        double beta0 = kp1*sin(thetai*deg)*sin(phi+zeta);
        double gamma00 = kp1*cos(thetai*deg);

        double muk0 = mu*k0;
        double muk0k0 = mu*sqr(k0);
        double seczetaovermuk0 = seczeta/muk0;


        Vr.allocate(M1,M2);
        Vt.allocate(M1,M2);
        R.allocate(M1,M2);
        r.allocate(M1,M2);
        T.allocate(M1,M2);
        t.allocate(M1,M2);

        DFARRAY alpha(M1,1);
        DFARRAY beta(M2,1);
        CFARRAY gamma_p1(M1,M2);
        CFARRAY gamma_m1(M1,M2);
        CFARRAY F(M1,M2,2,M1,M2,2);
        CFARRAY FG(M1,M2,2,M1,M2,2);
        CFARRAY G(M1,M2,2,M1,M2,2);

        CFARRAY W11(M1,M2,2,MM2);
        CFARRAY W21(M1,M2,2,MM2);
        CFARRAY W11inv(MM2,MM2);
        CFARRAY W21inv(MM2,MM2);
        CFARRAY W11old(MM2,MM2);
        CFARRAY W21old(MM2,MM2);

        CFARRAY t11(MM2,MM2);
        CFARRAY t21(MM2,MM2);

        CFARRAY Rud(MM2,MM2);
        CFARRAY Tdd(MM2,MM2);

        CFARRAY phip(MM2,1);
        CFARRAY t22_plus_t21Omega(MM2,MM2);
        CFARRAY t12_plus_t11Omega(MM2,MM2);
        CFARRAY gamma(MM2,1);
        CFARRAY E(M1,M2,2,MM2);
        CFARRAY H(M1,M2,2,MM2);

        CFARRAY work(3*MM2,1);
        IFARRAY pivot(MM2,1);
        CFARRAY det(2,1);
        int info;

        message("Approximate amount of memory allocated to RCW simulation: "
                + to_string((M1+M2+M1*M2*2+levels+M1*M2*2*M1*M2*2*16)*sizeof(COMPLEX)) + " bytes\n");



        for (i=1; i<=M1; ++i) {
            alpha(i) = alpha0 + (i-order1-1)*K1;
        }

        for (i=1; i<=M2; ++i) {
            beta(i) = beta0 + (i-order2-1)*K2;
        }



        for (i=1; i<=M1; ++i) {
            for (j=1; j<=M2; ++j) {

                COMPLEX kxy2=(sqr(alpha(i))+sqr(beta(j))-2*alpha(i)*beta(j)*sinzeta)*sqr(seczeta);
                COMPLEX g1 = sqrt(sqr(kp1)-kxy2);
                COMPLEX g2 = sqrt(sqr(km1)-kxy2);

                gamma_p1(i,j) = (real(g1)+imag(g1)>0) ? g1 : -g1;
                gamma_m1(i,j) = (real(g2)+imag(g2)>0) ? g2 : -g2;
            }
        }


        for (i=1; i<=MM2*MM2; ++i) Rud(i)=0.;
        for (i=1; i<=MM2*MM2; ++i) Tdd(i)=0.;
        for (i=1; i<=MM2; ++i) Tdd(i,i)=1.;




        CFARRAY _W11old=W11old;
        _W11old.array(M1,M2,2,M1,M2,2);
        CFARRAY _W21old=W21old;
        _W21old.array(M1,M2,2,M1,M2,2);

        COMPLEX em1muk0k0 = em1*mu*k0*k0;
        COMPLEX em1muk0k0sinzeta = em1muk0k0*sinzeta;

        for (i=1; i<=MM2*MM2; ++i) {
            _W11old(i)=0.;
            _W21old(i)=0.;
        }
        for (i=1; i<=M1; ++i) {
            for (j=1; j<=M2; ++j) {

                _W11old(i,j,1,i,j,1) = 1.;
                _W11old(i,j,1,i,j,2) = 0.;
                _W11old(i,j,2,i,j,1) = 0.;
                _W11old(i,j,2,i,j,2) = 1.;


                _W21old(i,j,1,i,j,1) = (em1muk0k0sinzeta-alpha(i)*beta(j))*seczetaovermuk0/gamma_m1(i,j);
                _W21old(i,j,1,i,j,2) = (alpha(i)*alpha(i)-em1muk0k0)*seczetaovermuk0/gamma_m1(i,j);
                _W21old(i,j,2,i,j,1) = (em1muk0k0-beta(j)*beta(j))*seczetaovermuk0/gamma_m1(i,j);
                _W21old(i,j,2,i,j,2) = (alpha(i)*beta(j)-em1muk0k0sinzeta)*seczetaovermuk0/gamma_m1(i,j);

            }
        }



        for (i=1; i<=MM2; ++i) {
            phip(i) = 1.;
        }




        ofstream eigenfile;
        if (write_eigenvalues.size()!=0) {
            eigenfile.open(write_eigenvalues.c_str());
            if (!eigenfile) error("Cannot open file " + write_eigenvalues);
        }




        forloopint levelloop(1,levels,(type==0||type==1) ? +1 : -1);
        for (forloopint::iterator level=levelloop.begin(); level!=levelloop.end(); ++level) {


            bool homogeneous = true;
            for (i=1; i<=M1; ++i) {
                for (j=1; j<=M2; ++j) {
                    for (k=1; k<=M1; ++k) {
                        for (l=1; l<=M2; ++l) {
                            if (i!=k||j!=l) {
                                if (abs(E0(i,j,k,l,level)) > 1E-10) {
                                    homogeneous = false;
                                }
                            }
                        }
                    }
                }
            }


            if (homogeneous) {

                COMPLEX elevel1 = E11(1,1,1,1,level);
				COMPLEX elevel2 = E12(1, 1, 1, 1, level);
				COMPLEX elevel3 = E0(1, 1, 1, 1, level);
				COMPLEX klevel1 = k0*sqrt(elevel1*mu);
				COMPLEX klevel2 = k0*sqrt(elevel2*mu);
				COMPLEX klevel3 = k0*sqrt(elevel3*mu);

				CFARRAY _gamma1 = gamma;
				CFARRAY _gamma2 = gamma;
				CFARRAY _gamma3 = gamma;
				_gamma1.array(M1,M2,2);
				_gamma2.array(M1, M2, 2);
				_gamma3.array(M1, M2, 2);

                k=0;
                for (i=1; i<=M1; ++i) {
                    for (j=1; j<=M2; ++j) {

                        COMPLEX kxy2 = (sqr(alpha(i))+sqr(beta(j))-2*alpha(i)*beta(j)*sinzeta)*sqr(seczeta);
						COMPLEX g1 = sqrt(sqr(klevel1) - kxy2);
						COMPLEX g2 = sqrt(sqr(klevel2) - kxy2);
						COMPLEX g3 = sqrt(sqr(klevel3) - kxy2);

                        g1 = (real(g1) + imag(g1)>0) ? g1 : -g1;
						g2 = (real(g2) + imag(g2)>0) ? g2 : -g2;
						g3 = (real(g3) + imag(g3)>0) ? g3 : -g3;
						_gamma1(i, j, 1) = g1;
						_gamma2(i, j, 1) = g2;
						_gamma3(i, j, 1) = g3;
						_gamma1(i, j, 2) = _gamma1(i, j, 1);
						_gamma2(i, j, 2) = _gamma2(i, j, 1);
						_gamma3(i, j, 2) = _gamma3(i, j, 1);
						++k;
                    }
                }


                if (write_eigenvalues.size()!=0) {
                    for (i=1; i<=M1; ++i) {
                        for (j=1; j<=M2; ++j) {
                            COMPLEX g1 = _gamma1(i,j,1);
                            eigenfile << level << tab << k << tab << real(g1) << tab << imag(g1) << endl;
                            eigenfile << level << tab << k << tab << real(g1) << tab << imag(g1) << endl;
                        }
                    }
                    eigenfile << endl;
                }

                CFARRAY _W11inv=W11inv;
                _W11inv.array(M1,M2,2,M1,M2,2);
                CFARRAY _W21inv=W21inv;
                _W21inv.array(M1,M2,2,M1,M2,2);
                CFARRAY _W11=W11;
                _W11.array(M1,M2,2,M1,M2,2);
                CFARRAY _W21=W21;
                _W21.array(M1,M2,2,M1,M2,2);


                COMPLEX elevel1muk0k0 = elevel1*mu*k0*k0;
                COMPLEX elevel1muk0k0sinzeta = elevel1muk0k0*sinzeta;
				COMPLEX elevel2muk0k0 = elevel2*mu*k0*k0;
				COMPLEX elevel2muk0k0sinzeta = elevel2muk0k0*sinzeta;
				COMPLEX elevel3muk0k0 = elevel3*mu*k0*k0;
				COMPLEX elevel3muk0k0sinzeta = elevel3muk0k0*sinzeta;

                for (i=1; i<=MM2*MM2; ++i) _W11inv(i)=0.;
                for (i=1; i<=MM2*MM2; ++i) _W21inv(i)=0.;
                for (i=1; i<=MM2*MM2; ++i) _W11(i)=0.;
                for (i=1; i<=MM2*MM2; ++i) _W21(i)=0.;
                for (i=1; i<=M1; ++i) {
                    for (j=1; j<=M2; ++j) {

                        _W11inv(i,j,1,i,j,1) = 1.;
                        _W11inv(i,j,1,i,j,2) = 0.;
                        _W11inv(i,j,2,i,j,1) = 0.;
                        _W11inv(i,j,2,i,j,2) = 1.;
                        _W11(i,j,1,i,j,1) = 1.;
                        _W11(i,j,1,i,j,2) = 0.;
                        _W11(i,j,2,i,j,1) = 0.;
                        _W11(i,j,2,i,j,2) = 1.;



                        COMPLEX a = (elevel1muk0k0sinzeta-alpha(i)*beta(j))*seczetaovermuk0/_gamma3(i,j);
                        COMPLEX b = (alpha(i)*alpha(i)-elevel2muk0k0)*seczetaovermuk0/_gamma3(i,j);
                        COMPLEX c = (elevel1muk0k0-beta(j)*beta(j))*seczetaovermuk0/_gamma3(i,j);

                        COMPLEX d = (alpha(i)*beta(j)-elevel2muk0k0sinzeta)*seczetaovermuk0/_gamma3(i,j);
                        COMPLEX det = a*d-b*c;
                        _W21inv(i,j,1,i,j,1) =  d/det;
                        _W21inv(i,j,1,i,j,2) = -b/det;
                        _W21inv(i,j,2,i,j,1) = -c/det;
                        _W21inv(i,j,2,i,j,2) =  a/det;
                        _W21(i,j,1,i,j,1) = a;
                        _W21(i,j,1,i,j,2) = b;
                        _W21(i,j,2,i,j,1) = c;
                        _W21(i,j,2,i,j,2) = d;

                    }
                }

            } else {


                for (i=1; i<=M1; ++i) {
                    for (j=1; j<=M2; ++j) {
                        for (k=1; k<=M1; ++k) {
                            for (l=1; l<=M2; ++l) {
                                COMPLEX _E0 = E0(i,j,k,l,level);
                                COMPLEX _E11 = E11(i,j,k,l,level);
								COMPLEX _E12 = E12(i,j,k,l,level);
								COMPLEX _E2 = E2(i,j,k,l,level);
                                COMPLEX _E3 = E3(i,j,k,l,level);
								COMPLEX _MU0 = MU0(i, j, k, l, level);
								COMPLEX _MU11 = MU11(i, j, k, l, level);
								COMPLEX _MU12 = MU12(i, j, k, l, level);
								COMPLEX _MU2 = MU2(i, j, k, l, level);
								COMPLEX _MU3 = MU3(i, j, k, l, level);

								COMPLEX sinzeta_E11 = sinzeta*_E11;
								COMPLEX sinzeta_E12 = sinzeta*_E12;
								COMPLEX sqrsinzeta_E11 = sinzeta*sinzeta_E11;
								COMPLEX sqrsinzeta_E12 = sinzeta*sinzeta_E12;
								COMPLEX sinzeta_MU11 = sinzeta*_MU11;
								COMPLEX sinzeta_MU12 = sinzeta*_MU12;
								COMPLEX sqrsinzeta_MU11 = sinzeta*sinzeta_MU11;
								COMPLEX sqrsinzeta_MU12 = sinzeta*sinzeta_MU12;
								F(i,j,1,k,l,1) =  alpha(i)*_E0*beta(l) - muk0k0*sinzeta_MU11;
                                F(i,j,1,k,l,2) = -alpha(i)*_E0*alpha(k) + muk0k0*(sqrcoszeta*_MU3+sqrsinzeta_MU12);
                                F(i,j,2,k,l,1) =  beta(j)* _E0*beta(l) - muk0k0*(sqrcoszeta*_MU2 + sqrsinzeta_MU11);
                                F(i,j,2,k,l,2) = -beta(j)* _E0*alpha(k) + muk0k0*sinzeta*_MU12;
								G(i,j,1,k,l,1) = -alpha(i)*_MU0*beta(l) + muk0k0*sinzeta_E11;
                                G(i,j,1,k,l,2) =  alpha(i)*_MU0*alpha(k) -muk0k0*(sqrcoszeta*_E3+sqrsinzeta_E12) ;
                                G(i,j,2,k,l,1) = -beta(j)* _MU0*beta(l)+ muk0k0*(sqrcoszeta*_E2+sqrsinzeta_E11);
								G(i,j,2,k,l,2) =  beta(j)*_MU0*alpha(k) -muk0k0*sinzeta_E12 ;
                            }
                        }
                    }
                }



                CMLIB::ZGEMM("N","N",MM2,MM2,MM2,1.,F,MM2,G,MM2,0.,FG,MM2);


                CMLIB::CGEEV((double*)&(FG[0]),MM2,MM2,(double*)&(gamma[0]),(double*)&(E[0]),MM2,(double*)&(work[0]),1,info);
                if (info!=0) error("Eigenvalue solver failed to converge");


                double muk0k0sqrcoszeta = muk0k0*sqr(coszeta);
                k=0;
                for (q=1; q<=MM2; ++q) {
                    COMPLEX g = sqrt(gamma(q)/muk0k0sqrcoszeta);
                    if (real(g) + imag(g)<0) g = -g;
                    gamma(q) = g;
                    eigenfile << level << tab << k << tab << real(g) << tab << imag(g) << endl;
                    eigenfile << level << tab << k << tab << real(g) << tab << imag(g) << endl;
                    ++k;
                }


                for (q=1; q<=MM2; ++q) {
                    COMPLEX seczetaovermuk0gamma = seczeta/muk0/gamma(q);
                    CMLIB::ZCOPY(MM2,&E(1,1,1,q),1,&H(1,1,1,q),1);
                    CMLIB::ZGEMV("N",MM2,MM2,seczetaovermuk0gamma,G,MM2,&E(1,1,1,q),1,0.,&H(1,1,1,q),1);
                }





                CMLIB::ZCOPY(MM2*MM2,E,1,W11,1);
                CMLIB::ZCOPY(MM2*MM2,H,1,W21,1);


                CMLIB::ZCOPY(MM2*MM2,W11,1,W11inv,1);
                CMLIB::ZCOPY(MM2*MM2,W21,1,W21inv,1);


                CMLIB::CGEFA(W11inv,MM2,MM2,pivot,info);
                if (info!=0) error("Singular matrix in LUdecompose");
                CMLIB::CGEDI(W11inv,MM2,MM2,pivot,det,work,1);

                CMLIB::CGEFA(W21inv,MM2,MM2,pivot,info);
                if (info!=0) error("Singular matrix in LUdecompose");
                CMLIB::CGEDI(W21inv,MM2,MM2,pivot,det,work,1);

            }





            for (i=1; i<=MM2; ++i) {
                for (j=1; j<=MM2; ++j) {
                    t11(i,j) = 0;
                    t21(i,j) = 0;
                    Rud(i,j) = phip(i)*Rud(i,j)*phip(j);
                    Tdd(i,j) = Tdd(i,j)*phip(j);
                    for (k=1; k<=MM2; ++k) {
                        COMPLEX AA = W11inv(i,k)*W11old(k,j)/2.;
                        COMPLEX BB = W21inv(i,k)*W21old(k,j)/2.;
                        t11(i,j) += AA+BB;
                        t21(i,j) += AA-BB;
                    }
                }
            }


            CMLIB::ZCOPY(MM2*MM2,W11,1,W11old,1);
            CMLIB::ZCOPY(MM2*MM2,W21,1,W21old,1);




            CMLIB::ZCOPY(MM2*MM2,t11,1,t22_plus_t21Omega,1);
            CMLIB::ZCOPY(MM2*MM2,t21,1,t12_plus_t11Omega,1);

            CMLIB::ZGEMM("N","N",MM2,MM2,MM2,1.,t21,MM2,Rud,MM2,1.,t22_plus_t21Omega,MM2);
            CMLIB::ZGEMM("N","N",MM2,MM2,MM2,1.,t11,MM2,Rud,MM2,1.,t12_plus_t11Omega,MM2);

            for (i=1; i<=MM2; ++i) {
                for (j=i+1; j<=MM2; ++j) {
                    SWAP(t22_plus_t21Omega(i,j),t22_plus_t21Omega(j,i));
                    SWAP(t12_plus_t11Omega(i,j),t12_plus_t11Omega(j,i));
                    SWAP(Tdd(i,j),Tdd(j,i));
                }
            }


            CMLIB::CGEFA(t22_plus_t21Omega,MM2,MM2,pivot,info);
            if (info!=0) error("Singular matrix in LUdecompose");

            CMLIB::ZCOPY(MM2*MM2,t12_plus_t11Omega,1,Rud,1);

            for (i=1; i<=MM2; ++i) {
                CMLIB::CGESL(t22_plus_t21Omega,MM2,MM2,pivot,Rud(1,i),0);
                CMLIB::CGESL(t22_plus_t21Omega,MM2,MM2,pivot,Tdd(1,i),0);
            }

            for (i=1; i<=MM2; ++i) {
                for (j=i+1; j<=MM2; ++j) {
                    SWAP(Rud(i,j),Rud(j,i));
                    SWAP(Tdd(i,j),Tdd(j,i));
                }
            }




            for (i=1; i<=MM2; ++i) {
                phip(i) = exp(cI*gamma(i)*thick(level));
            }

        }

        CFARRAY _W11inv=W11inv;
        _W11inv.array(M1,M2,2,M1,M2,2);
        CFARRAY _W21inv=W21inv;
        _W21inv.array(M1,M2,2,M1,M2,2);


        COMPLEX ep1muk0k0 = ep1*mu*k0*k0;
        COMPLEX ep1muk0k0sinzeta = ep1muk0k0*sinzeta;

        for (i=1; i<=MM2*MM2; ++i) _W11inv(i)=0.;
        for (i=1; i<=MM2*MM2; ++i) _W21inv(i)=0.;
        for (i=1; i<=M1; ++i) {
            for (j=1; j<=M2; ++j) {

                _W11inv(i,j,1,i,j,1) = 1.;
                _W11inv(i,j,1,i,j,2) = 0.;
                _W11inv(i,j,2,i,j,1) = 0.;
                _W11inv(i,j,2,i,j,2) = 1.;


                COMPLEX a = (ep1muk0k0sinzeta-alpha(i)*beta(j))*seczetaovermuk0/gamma_p1(i,j);
                COMPLEX b = (alpha(i)*alpha(i)-ep1muk0k0)*seczetaovermuk0/gamma_p1(i,j);
                COMPLEX c = (ep1muk0k0-beta(j)*beta(j))*seczetaovermuk0/gamma_p1(i,j);
                COMPLEX d = (alpha(i)*beta(j)-ep1muk0k0sinzeta)*seczetaovermuk0/gamma_p1(i,j);
                COMPLEX det = a*d-b*c;
                _W21inv(i,j,1,i,j,1) =  d/det;
                _W21inv(i,j,1,i,j,2) = -b/det;
                _W21inv(i,j,2,i,j,1) = -c/det;
                _W21inv(i,j,2,i,j,2) =  a/det;
            }
        }






        for (i=1; i<=MM2; ++i) {
            for (j=1; j<=MM2; ++j) {
                t11(i,j) = 0;
                t21(i,j) = 0;
                Rud(i,j) = phip(i)*Rud(i,j)*phip(j);
                Tdd(i,j) = Tdd(i,j)*phip(j);
                for (k=1; k<=MM2; ++k) {
                    COMPLEX AA = W11inv(i,k)*W11old(k,j)/2.;
                    COMPLEX BB = W21inv(i,k)*W21old(k,j)/2.;
                    t11(i,j) += AA+BB;
                    t21(i,j) += AA-BB;
                }
            }
        }




        CMLIB::ZCOPY(MM2*MM2,t11,1,t22_plus_t21Omega,1);
        CMLIB::ZCOPY(MM2*MM2,t21,1,t12_plus_t11Omega,1);

        CMLIB::ZGEMM("N","N",MM2,MM2,MM2,1.,t21,MM2,Rud,MM2,1.,t22_plus_t21Omega,MM2);
        CMLIB::ZGEMM("N","N",MM2,MM2,MM2,1.,t11,MM2,Rud,MM2,1.,t12_plus_t11Omega,MM2);

        for (i=1; i<=MM2; ++i) {
            for (j=i+1; j<=MM2; ++j) {
                SWAP(t22_plus_t21Omega(i,j),t22_plus_t21Omega(j,i));
                SWAP(t12_plus_t11Omega(i,j),t12_plus_t11Omega(j,i));
                SWAP(Tdd(i,j),Tdd(j,i));
            }
        }


        CMLIB::CGEFA(t22_plus_t21Omega,MM2,MM2,pivot,info);
        if (info!=0) error("Singular matrix in LUdecompose");

        CMLIB::ZCOPY(MM2*MM2,t12_plus_t11Omega,1,Rud,1);

        for (i=1; i<=MM2; ++i) {
            CMLIB::CGESL(t22_plus_t21Omega,MM2,MM2,pivot,Rud(1,i),0);
            CMLIB::CGESL(t22_plus_t21Omega,MM2,MM2,pivot,Tdd(1,i),0);
        }

        for (i=1; i<=MM2; ++i) {
            for (j=i+1; j<=MM2; ++j) {
                SWAP(Rud(i,j),Rud(j,i));
                SWAP(Tdd(i,j),Tdd(j,i));
            }
        }


        CFARRAY _Rud(Rud);
        _Rud.array(M1,M2,2,M1,M2,2);

        CFARRAY _Tdd(Tdd);
        _Tdd.array(M1,M2,2,M1,M2,2);











        CVector bsub1(1,0,0);
        CVector bsub2(sinzeta,coszeta,0);
        CVector bsub3(0,0,1);

        CVector bsup1(1,-tan(zeta),0);
        CVector bsup2(0,seczeta,0);
        CVector bsup3(0,0,1);



        CVector Sin(-sinphi,cosphi,0);
        CVector Pin(cos(thetai*deg)*cosphi,cos(thetai*deg)*sinphi,sin(thetai*deg));
        CVector Kin(kp1*sin(thetai*deg)*cosphi,kp1*sin(thetai*deg)*sinphi,kp1*cos(thetai*deg));

        CVector zhat(0,0,1);

        COMPLEX Ei1s = Sin*bsub1;
        COMPLEX Ei2s = Sin*bsub2;
        COMPLEX Ei1p = Pin*bsub1;
        COMPLEX Ei2p = Pin*bsub2;

        for (i=1; i<=M1; ++i) {
            for (j=1; j<=M2; ++j) {


                if (1) {

                    CVector a(alpha(i),beta(j),gamma_p1(i,j));


                    CVector Kout(a.x,-a.x*tanzeta+a.y*seczeta,a.z);
                    COMPLEX rho = sqrt(sqr(Kout.x)+sqr(Kout.y));
                    CVector Sout = norm(rho)!=0 ? CVector(-Kout.y/rho,Kout.x/rho,0) : CVector(-sinphi,cosphi,0);
                    CVector Pout = cross(Kout,Sout)/kp1;


                    Vr(i,j) = CVector(Kout.x*cosphi+Kout.y*sinphi,-Kout.x*sinphi+Kout.y*cosphi,Kout.z);

                    COMPLEX Er1s = bsup1*Sout;
                    COMPLEX Er2s = bsup2*Sout;
                    COMPLEX Er1p = bsup1*Pout;
                    COMPLEX Er2p = bsup2*Pout;

                    int o1 = order1+1;
                    int o2 = order2+1;
                    COMPLEX R11 = _Rud(i,j,1,o1,o2,1);
                    COMPLEX R12 = _Rud(i,j,1,o1,o2,2);
                    COMPLEX R21 = _Rud(i,j,2,o1,o2,1);
                    COMPLEX R22 = _Rud(i,j,2,o1,o2,2);

                    COMPLEX A1s = R11*Ei1s + R12*Ei2s;
                    COMPLEX A2s = R21*Ei1s + R22*Ei2s;
                    COMPLEX A1p = R11*Ei1p + R12*Ei2p;
                    COMPLEX A2p = R21*Ei1p + R22*Ei2p;

                    COMPLEX ct = Kout.z/kp1;
                    COMPLEX Rss = A1s*Er1s + A2s*Er2s;
                    COMPLEX Rsp = A1p*Er1s + A2p*Er2s;
                    COMPLEX Rps = (A1s*Er1p + A2s*Er2p)/sqr(ct);
                    COMPLEX Rpp = (A1p*Er1p + A2p*Er2p)/sqr(ct);

                    JonesMatrix jones;
                    if (type==0 || type==1) {
                        jones.SS() =  Rss;
                        jones.SP() =  Rps;
                        jones.PS() =  Rsp;
                        jones.PP() =  Rpp;
                    } else {
                        jones.SS() =  Rss;
                        jones.SP() = -Rps;
                        jones.PS() = -Rsp;
                        jones.PP() =  Rpp;
                    }
                    r(i,j) = jones;
                    R(i,j) = MuellerMatrix(jones)/fabs(real(gamma00/gamma_p1(i,j)));
                }

                if (1) {

                    CVector b(alpha(i),beta(j),-gamma_m1(i,j));


                    CVector Kout(b.x,-b.x*tanzeta+b.y*seczeta,b.z);
                    COMPLEX rho = sqrt(sqr(Kout.x)+sqr(Kout.y));
                    CVector Sout = norm(rho)!=0 ? CVector(-Kout.y/rho,Kout.x/rho,0) : CVector(-sinphi,cosphi,0);
                    CVector Pout = cross(Kout,Sout)/km1;


                    Vt(i,j) = CVector(Kout.x*cosphi+Kout.y*sinphi,-Kout.x*sinphi+Kout.y*cosphi,Kout.z);

                    COMPLEX Er1s = bsup1*Sout;
                    COMPLEX Er2s = bsup2*Sout;
                    COMPLEX Er1p = bsup1*Pout;
                    COMPLEX Er2p = bsup2*Pout;

                    int o1 = order1+1;
                    int o2 = order2+1;
                    COMPLEX T11 = _Tdd(i,j,1,o1,o2,1);
                    COMPLEX T12 = _Tdd(i,j,1,o1,o2,2);
                    COMPLEX T21 = _Tdd(i,j,2,o1,o2,1);
                    COMPLEX T22 = _Tdd(i,j,2,o1,o2,2);

                    COMPLEX A1s = T11*Ei1s + T12*Ei2s;
                    COMPLEX A2s = T21*Ei1s + T22*Ei2s;
                    COMPLEX A1p = T11*Ei1p + T12*Ei2p;
                    COMPLEX A2p = T21*Ei1p + T22*Ei2p;

                    COMPLEX ct = Kout.z/km1;
                    COMPLEX Tss = A1s*Er1s + A2s*Er2s;
                    COMPLEX Tsp = A1p*Er1s + A2p*Er2s;
                    COMPLEX Tps = (A1s*Er1p + A2s*Er2p)/sqr(ct);
                    COMPLEX Tpp = (A1p*Er1p + A2p*Er2p)/sqr(ct);

                    JonesMatrix jones;
                    if (type==0 || type==1) {
                        jones.SS() =  Tss;
                        jones.SP() =  Tps;
                        jones.PS() =  Tsp;
                        jones.PP() =  Tpp;
                    } else {
                        jones.SS() =  Tss;
                        jones.SP() = -Tps;
                        jones.PS() = -Tsp;
                        jones.PP() =  Tpp;
                    }
                    t(i,j) = jones;
                    T(i,j) = MuellerMatrix(jones)/fabs(real(gamma00/gamma_m1(i,j)));
                }
            }
        }
        old_type = type;
    }

    Vector CrossRCW_Model::GetDirection(int i,int j)
    {
        SETUP();

        CVector _V = GetPropagationVector(i,j);

        if (imag(_V.z)!=0.) return Vector(0,0,0);

        Vector result(real(_V.x),real(_V.y),real(_V.z));
        return result/Norm(result);
    }

    CVector CrossRCW_Model::GetPropagationVector(int i,int j)
    {
        SETUP();

        i = i+order1+1;
        j = j+order2+1;

        CVector result = ((type&0x01) ? Vt(i,j) : Vr(i,j));
		if (type == 2 || type == 3) result.z = -result.z;
		return result;
    }

    JonesMatrix CrossRCW_Model::GetAmplitude(int i,int j)
    {
        SETUP();

        i = i+order1+1;
        j = j+order2+1;

        return (type&0x01) ? t(i,j) : r(i,j);
    }

    MuellerMatrix CrossRCW_Model::GetIntensity(int i,int j)
    {
        SETUP();

        i = i+order1+1;
        j = j+order2+1;

        return (type&0x01) ? T(i,j) : R(i,j);
    }

	StokesVector CrossRCW_Model::GetAbsorption()
	{
		SETUP();
		StokesVector result(0,0,0,0);
		for (int i=1;i<=2*order1+1;++i) {
			for (int j=1;j<=2*order2+1;++j) {
				result += (T(i,j) + R(i,j)).transpose()*StokesVector(1,0,0,0);
			}
		}
		return StokesVector(1,0,0,0)-result;
	}

    CVector CrossRCW_Model::GetEField(const JonesVector& inpol, const Vector& pos, bool incident)
    {
        SETUP();

        CVector result(0.,0.,0.);

        if (type==0 && pos.z<0) return result;
        if (type==1 && pos.z>-totalthickness) return result;
        if (type==2 && pos.z>-totalthickness) return result;
        if (type==3 && pos.z<0) return result;

        Vector R = pos;
        if (type==1) R.z += totalthickness;
        if (type==2) R.z += totalthickness;

        if (type==0 && incident) {
            Vector kin(nI*k0*sin(thetai*deg),0.,-nI*k0*cos(thetai*deg));
            Vector s(0.,1.,0.);
            Vector p = cross(kin,s)/(k0*nI);
            COMPLEX phase = exp(cI*COMPLEX(kin*R));
            result = (inpol.S()*CVector(s)+inpol.P()*CVector(p))*phase;
        }

        if (type==2 && incident) {
            CVector kin(conj(nII)*k0*sin(thetai*deg),0.,conj(nII)*k0*cos(thetai*deg));
            CVector s(0.,1.,0.);
            CVector p = cross(kin,s)/(k0*conj(nII));
            COMPLEX phase = exp(cI*COMPLEX(kin*(CVector)R));
            result = (inpol.S()*CVector(s)+inpol.P()*CVector(p))*phase;
        }

        for (int i=-order1; i<=order1; ++i) {
            for (int j=-order2; j<=order2; ++j) {
                JonesVector out = GetAmplitude(i,j)*inpol;
                CVector kvector = GetPropagationVector(i,j);
                double x = real(kvector.x);
                double y = real(kvector.y);
                double rho = sqrt(x*x+y*y);
                COMPLEX k = k0 * (type==0 || type==3 ? nI : conj(nII));

				CVector s = CVector(-y / rho, x / rho, 0);
                CVector p = cross(kvector,s)/k;
                COMPLEX phase = exp(cI*(kvector*(CVector)R));
                result += (out.S()*s+out.P()*p)*phase;
            }
        }
        return result;
    }

    CVector CrossRCW_Model::GetHField(const JonesVector& inpol, const Vector& pos, bool incident)
    {
        SETUP();

        CVector result(0.,0.,0.);

        if (type==0 && pos.z<0) return result;
        if (type==1 && pos.z>-totalthickness) return result;
        if (type==2 && pos.z>-totalthickness) return result;
        if (type==3 && pos.z<0) return result;

        Vector R = pos;
        if (type==1) R.z += totalthickness;
        if (type==2) R.z += totalthickness;

        if (type==0 && incident) {
            Vector kin(nI*k0*sin(thetai*deg),0.,-nI*k0*cos(thetai*deg));
            Vector p(0.,-1,0.);
            Vector s = cross(p,kin)/(k0*nI);
            COMPLEX phase = nI*exp(cI*COMPLEX(kin*R));
            result = (inpol.S()*CVector(s)+inpol.P()*CVector(p))*phase;
        }

        if (type==2 && incident) {
            CVector kin(conj(nII)*k0*sin(thetai*deg),0.,conj(nII)*k0*cos(thetai*deg));
            CVector p(0.,-1.,0.);
            CVector s = cross(p,kin)/(k0*conj(nII));
            COMPLEX phase = conj(nII)*exp(cI*COMPLEX(kin*(CVector)R));
            result = (inpol.S()*CVector(s)+inpol.P()*CVector(p))*phase;
        }

        for (int i=-order1; i<=order1; ++i) {
            for (int j=-order2; j<=order2; ++j) {
                JonesVector out = GetAmplitude(i,j)*inpol;
                CVector kvector = GetPropagationVector(i,j);
                double x = real(kvector.x);
                double y = real(kvector.y);
                COMPLEX N = (type==0 || type==3 ? nI : conj(nII));
                double rho = sqrt(x*x+y*y);
                COMPLEX k = k0 * N;

				CVector p = CVector(y / rho, -x / rho, 0);
                CVector s = cross(p,kvector)/k;
                COMPLEX phase = N*exp(cI*(kvector*(CVector)R));
                result += (out.S()*s+out.P()*p)*phase;
            }
        }
        return result;
    }

    CVector CrossRCW_Model::GetBField(const JonesVector& inpol, const Vector& pos, bool incident)
    {
        return GetHField(inpol,pos,incident);
    }

    CVector CrossRCW_Model::GetDField(const JonesVector& inpol, const Vector& pos, bool incident)
    {
        return GetEField(inpol,pos,incident) * (type==0 || type==3 ? sqr(nI) : sqr(conj(nII)));
    }


    void CrossRCW_Model::set_parameter_base(
        const STRING& parameter,
        const STRING& value
    )
    {
        if (parameter=="WriteEigenvalues") {
            write_eigenvalues=value;
            set_recalc(-1);
        } else {
            Model::set_parameter_base(parameter,value);
        }
    }

    void CrossRCW_BRDF_Model::setup()
    {
        BRDF_Model::setup();

        if (RCW.get_lambda()!=lambda) RCW.set_lambda(lambda);
        if (RCW.get_thetai()!=thetai/deg) RCW.set_thetai(thetai/deg);
        if (RCW.get_rotation()!=rotation/deg) RCW.set_rotation(rotation/deg);
        if (RCW.get_order1()!=order1) RCW.set_order1(order1);
        if (RCW.get_order2()!=order2) RCW.set_order2(order2);
        if (RCW.get_type()!=type) RCW.set_type(type);
        if (grating->get_medium_t().index(lambda)!=substrate.index(lambda)) error("grating.medium_t!=substrate");
        if (grating->get_medium_i().index(lambda)!=COMPLEX(1,0)) error("grating.medium_i!=vacuum");
        if (grating->get_lambda()!=lambda) grating->set_lambda(lambda);
        if (grating->get_order1()!=order1) grating->set_order1(order1);
        if (grating->get_order2()!=order2) grating->set_order2(order2);
        if (grating->get_recalc()) {
            grating->SETUP();
            RCW.set_grating(grating);
        }

        M1 = 2*order1+1;
        M2 = 2*order2+1;

        V.allocate(M1,M2);
        R.allocate(M1,M2);

		for (int i=1; i<=M1; ++i) {
            for (int j=1; j<=M2; ++j) {
                V(i,j)=RCW.GetDirection(i-order1-1,j-order2-1);
                R(i,j)=RCW.GetIntensity(i-order1-1,j-order2-1);
			}
        }

        if (alpha==0.) error("alpha cannot be zero");
        cosalpha = cos(alpha);
        Omega = pi*sqr(alpha);
    }

    MuellerMatrix CrossRCW_BRDF_Model::mueller()
    {
        SETUP();
        Vector v = polar(1.,thetas,phis);
        MuellerMatrix result = MuellerZero();

        for (int i=1; i<=M1; ++i) {
            for (int j=1; j<=M2; ++j) {
                Vector _V = V(i,j);
				if (type == 1 || type == 2) _V.z = -_V.z;
                if (_V.z!=0) {
                    if (_V*v>=cosalpha) {
                        result += R(i,j)/(Omega*v.z);
                    }
                }
            }
        }
        return result;
    }

    DEFINE_MODEL(CrossRCW_Model,Model,"RCW solution for crossed gratings");
    DEFINE_PARAMETER(CrossRCW_Model,double,thetai,"Incident angle [deg]","0",0x01);
    DEFINE_PARAMETER(CrossRCW_Model,double,rotation,"Rotation angle of sample [deg]","0",0x01);
    DEFINE_PARAMETER(CrossRCW_Model,double,lambda,"Wavelength [um]","0.5",0x01);
    DEFINE_PARAMETER(CrossRCW_Model,int,type,"(0) for Forward/Reflection, (1) for Forward/Transmission, (2) for Backward/Reflection, or (3) for Backward/Transmission","0",0x02);
    DEFINE_PARAMETER(CrossRCW_Model,int,order1,"Fourier order in direction #1","10",0x01);
    DEFINE_PARAMETER(CrossRCW_Model,int,order2,"Fourier order in direction #2","10",0x01);
    DEFINE_PTRPARAMETER(CrossRCW_Model,CrossGrating_Ptr,grating,"Grating","OneD_CrossGrating",0x01);

    DEFINE_MODEL(CrossRCW_BRDF_Model,BRDF_Model,"RCW solution for crossed gratings");
    DEFINE_PARAMETER(CrossRCW_BRDF_Model,double,alpha,"Half angle of detector [rad]","0.01",0xFF);
    DEFINE_PARAMETER(CrossRCW_BRDF_Model,int,order1,"Fourier order in direction #1","10",0x01);
    DEFINE_PARAMETER(CrossRCW_BRDF_Model,int,order2,"Fourier order in direction #2","10",0x01);
    DEFINE_PTRPARAMETER(CrossRCW_BRDF_Model,CrossGrating_Ptr,grating,"Grating","OneD_CrossGrating",0xFF);

}

