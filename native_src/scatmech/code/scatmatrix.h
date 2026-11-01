












#ifndef SCATMECH_SCATMATRIX_H
#define SCATMECH_SCATMATRIX_H

#include "scatmech.h"

namespace SCATMECH {

    class ScatterTMatrix {
        public:
            ScatterTMatrix(int lmax=0);
            ScatterTMatrix(const ScatterTMatrix& matrix2);
            ScatterTMatrix& operator=(const ScatterTMatrix& matrix2);
            ~ScatterTMatrix();

            void resize(int lmax);
            void free();

            COMPLEX** operator[](int m) {
                return matrix[m];
            }
        private:
            COMPLEX ***matrix;
            int LMAX;
    };

}

#endif
