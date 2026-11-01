











#include "axipart.h"
#include <fstream>
#include <algorithm>

using namespace std;

namespace SCATMECH {

    using namespace BobVlieg_Supp;

    static void
    solve(vector<vector<COMPLEX> >& a,vector<vector<COMPLEX> >& b,int n)
    {




        vector<int> ls(n);



        for (int i=0; i<n; ++i) {



            COMPLEX aijmax = a[i][0];
            int jmax = 0;
            for (int j=1; j<n; ++j) {
                if(abs(a[i][j])>abs(aijmax)) {
                    aijmax = a[i][j];
                    jmax = j;
                }
            }



            for (int jj=0; jj<n; ++jj) {
                a[i][jj] /= aijmax;
                b[i][jj] /= aijmax;
            }





            for (int k=0; k<n; ++k) {
                if(k!=i) {
                    COMPLEX arat = -a[k][jmax];
                    for (int j=0; j<n; ++j) {
                        if(abs(a[i][j])>0.0) {
                            a[k][j] += arat*a[i][j];
                        }
                    }
                    a[k][jmax] = 0.0;
                    for (int jj=0; jj<n; ++jj) {
                        if(abs(b[i][jj])>0.0) {
                            b[k][jj] += arat*b[i][jj];
                        }
                    }
                }





                ls[jmax] = i;
            }
        }




        for (int ii=0; ii<n; ++ii) {
            int k = ii;








            do {
                k = ls[k];
            } while (k<ii);
            if (k>0) {
                for (int j=0; j<n; ++j) {
                    swap(b[ii][j],b[k][j]);
                }
            }
        }



        return;
    }

    void
    Axisymmetric_Shape::
    Get_TMatrix(ScatterTMatrix& T, int nrank, int mmax, double k, COMPLEX& index)
    {
































































        vector<double> pnmllg;
        vector<COMPLEX> hankel;
        vector<COMPLEX> bslcmp;

        int isize=2*(nrank+1);

        vector<vector<COMPLEX> > a(isize,vector<COMPLEX>(isize));
        vector<vector<COMPLEX> > b(isize,vector<COMPLEX>(isize));




        double pi2 = pi/2.0;

        int nranki = nrank+1;
        int nr2 = 2*nrank;





        COMPLEX cm = index;
        COMPLEX dcn = sqr(cm);





        int nm = (mmax<nrank) ? mmax : nrank;





        shape_vec shape2 = shape;

        int ntheta = shape2.size();

        for (int i=0; i<ntheta; ++i) {
            shape2[i].shape *= k;
            shape2[i].dshape *= k;
        }




        for (int im = 1; im<=nm+1; ++im) {
            ostringstream msg;
            msg << "Calculating T matrix: m = " << im-1 << ret;
            message(msg.str());




            int kmv = im-1;
            if (nm==1) kmv = 1;

            double cmv = kmv;
            double cm2 = sqr(cmv);
            double prodm = 1.0;

            double em;
            if(kmv==0) {
                em = 1.0;
            } else {
                em = 2.0;
                double quanm = cmv;
                for (int i=1; i<=kmv; ++i) {
                    quanm = quanm+1.0;
                    prodm = quanm*prodm/2.0;
                }
            }
            double qem = 2.0/em;
            double twm = 2.0*cmv;




            int ij = kmv-1;
            if (ij<0) ij = 0;
            int ijt = 2*ij;
            int ns = nrank-ij;
            int ns2 = 2*ns;



            for (int i=0; i<isize; ++i) {
                for (int j=0; j<isize; ++j) {
                    a[i][j] = b[i][j] = 0.;
                }
            }




            for (int ithta=0; ithta<ntheta; ++ithta) {

                double theta = shape2[ithta].abs;
                double sinth = sin(theta);
                double costh = cos(theta);





                pnmllg.resize(nranki+1);

                for (int n=0; n<=nranki; ++n) {
                    if (n>=kmv) {
                        pnmllg[n]=real(Legendre(n,kmv,costh))/sinth;
                    } else {
                        pnmllg[n]=0.;
                    }
                }
                pnmllg[nranki]=0.;





                double ckr = shape2[ithta].shape;
                double dckr = shape2[ithta].dshape;





                COMPLEX ckpr = cm*ckr;
                hankel.resize(nranki+1);
                bslcmp.resize(nranki+1);
                for (int ii=0; ii<nranki; ++ii) {
                    hankel[ii] = h(ii,ckr);
                    bslcmp[ii] = j(ii,ckpr);
                }
                hankel[nranki]=bslcmp[nranki]=0.;

                double d = dckr*sinth;
                double wtsin = shape2[ithta].wt*sinth;



                for (int irow=1; irow<=nrank; ++irow) {

                    if(irow>ij) {
                        int irow1 = irow+nrank;
                        double crow = irow;
                        double crowm = cmv+crow;
                        double crow1 = crow+1.0;
                        double br = real(hankel[irow-1])/real(hankel[irow+1-1]);
                        COMPLEX h = hankel[irow-1]/hankel[irow+1-1];



                        for (int icol=1; icol<=nrank; ++icol) {

                            if (icol>ij) {

                                int icol1 = icol+nrank;

                                double ccol = icol;
                                double ccolm = cmv+ccol;
                                double ccol1 = ccol+1.0;



                                double crij = crow+ccol;
                                double crssij = crow*ccol;
                                double cmcrco = cm2+qem*crssij*sqr(costh);
                                double pnr0c0 = pnmllg[irow-1]*pnmllg[icol-1];
                                double pnr0c1 = pnmllg[irow-1]*pnmllg[icol+1-1];
                                double pnr1c0 = pnmllg[irow+1-1]*pnmllg[icol-1];
                                double pnr1c1 = pnmllg[irow+1-1]*pnmllg[icol+1-1];
                                double b1a = crow*costh*pnr1c1-crowm*pnr0c1;
                                double b1b = ccol*costh*pnr1c1-ccolm*pnr1c0;
                                COMPLEX bk = cm*bslcmp[icol-1]/bslcmp[icol+1-1];
                                COMPLEX bbk = real(hankel[irow+1-1])*bslcmp[icol+1-1]*wtsin;
                                COMPLEX hbk = hankel[irow+1-1]*bslcmp[icol+1-1]*wtsin;


















                                if(kmv != 0) {
                                    double b1 = b1a+b1b;
                                    COMPLEX sa = pnr1c1*
                                                 (crow*crow1*bk+ccol*ccol1*h-crssij*(crij+2.0)/ckr)*d;
                                    sa = (ckr*(1.0+h*bk)-ccol*h-crow*bk+crssij/ckr)*b1*ckr+sa;

                                    a[icol-ij-1][irow1-ijt-1] += cmv*sa*hbk;
                                    COMPLEX sb = sa+(br-h)*(ccol*ccol1*pnr1c1*d+(bk*ckr-ccol)*b1*ckr);
                                    b[icol-ij-1][irow1-ijt-1] += +cmv*sb*bbk;



                                    COMPLEX c = (dcn-1.0)*b1*sqr(ckr);

                                    a[icol1-ijt-1][irow-ij-1] -= cmv*(sa+c)*hbk/cm;
                                    b[icol1-ijt-1][irow-ij-1] -= cmv*(sb+c)*bbk/cm;
                                }



                                double a12 = cmcrco*pnr1c1-qem*(crow*ccolm*costh*pnr1c0+
                                                                ccol*crowm*costh*pnr0c1-crowm*ccolm*pnr0c0);
                                double a22 = a12*sqr(ckr);
                                b1a = ccol*ccol1*b1a;
                                b1b = crow*crow1*b1b;
                                COMPLEX sa = (ckr*(bk-dcn*h)+dcn*crow-ccol)*a12*ckr+(b1a-dcn*b1b)*qem*d;

                                a[icol1-ijt-1][irow1-ijt-1] += sa*hbk/cm;
                                b[icol1-ijt-1][irow1-ijt-1] += (sa-(br-h)*dcn*a22)*bbk/cm;



                                sa = (ckr*(bk-h)+crow-ccol)*a12*ckr+(b1a-b1b)*qem*d;

                                a[icol-ij-1][irow-ij-1] += sa*hbk;
                                b[icol-ij-1][irow-ij-1] += (sa-(br-h)*a22)*bbk;

                            }
                        }
                    }
                }
            }




            solve(a,b,ns2);




            for (int ir=0; ir<ns2; ++ir) {
                for (int jr=0; jr<ns2; ++jr) {
                    a[ir][jr]=b[jr][ir];
                }
            }








            int lmax=nrank;
            int m = kmv;





            int mm = lmax+m;
            int mm2 = lmax-m;
            int beginl = (m==0) ? 1 : abs(m);
            for (int l=beginl; l<=lmax; ++l) {
                for (int f=0; f<=1; ++f) {
                    int ll = 2*(lmax-l)+f;
                    for (int l_=beginl; l_<=lmax; ++l_) {
                        for (int f_=0; f_<=1; ++f_) {
                            int ll_ = 2*(lmax-l_)+f_;

                            int nn = ((f==0) ? ns : 0 ) + (l-beginl);
                            int nn_ = ((f_==0) ? ns : 0 ) + (l_-beginl);

                            T[mm][ll_][ll]=-a[nn_][nn];
                            T[mm2][ll_][ll]= (f==f_) ? -a[nn_][nn] : a[nn_][nn];
                        }
                    }
                }
            }
        }
    }

}
