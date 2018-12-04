//
// Created by daniel.devries on 11/30/2018.
//

#include <cmath>
#include "xutils.h"

namespace xutils {

    /**
     * Solve general NxN system in N unknowns with arbitrary number of righthand sides.
     *
     * Asumes system is invertible...
     * ...if it isn't, a divide by zero will result.
     *
     * The matrix r represents the righthand side. If only one system is to be solved,
     * this should be a Nx1 matrix. If it is an NxM matrix, M systems of equations
     * with the same coefficient matrix will be solved.
     *
     * @param z     coefficient matrix (destroyed by the solution process)
     * @param r     right hand sides(s) (replaced by the solution vector(s))
     */
    template <int nsiz, int nrhs>
    void gauss(int nn, double z[nsiz][nsiz], double r[nsiz][nrhs]) {
        int np, np1, nx, n, l, k;
        double pivot, temp, ztmp;
        for (np = 0; np < nn - 1; np++) {
            np1 = np + 1;

            // find max pivot index nx
            nx = np;
            for (n = np1; n < nn; n++) {
                if (fabs(z[n][np]) - fabs(z[nx][np]) > 0) {
                    nx = n;
                }
            }

            pivot = 1.0 / z[nx][np];

            // switch pivots
            z[nx][np] = z[np][np];

            // switch rows & normalize pivot row
            for (l = np1; l < nn; l++) {
                temp = z[nx][l] * pivot;
                z[nx][l] = z[np][l];
                z[np][l] = temp;
            }

            for (l = 0; l < nrhs; l++) {
                temp = r[nx][l] * pivot;
                r[nx][l] = r[np][l];
                r[np][l] = temp;
            }

            // forward eliminate everything
            for (k = np1; k < nn; k++) {
                ztmp = z[k][np];

                for (l = np1; l < nn; l++) {
                    z[k][l] -= ztmp * z[np][l];
                }
                for (l = 0; l < nrhs; l++) {
                    r[k][l] -= ztmp * r[np][l];
                }
            }
        }

        // solve for last row
        for (l = 0; l < nrhs; l++) {
            r[nn-1][l] /= z[nn-1][nn-1];
        }

        // back substitute everything
        for (np = nn - 2; np >= 0; np--) {
            np1 = np + 1;
            for (l = 0; l < nrhs; l++) {
                for (k = np1; k < nn; k++) {
                    r[np][l] -= z[np][k] * r[k][l];
                }
            }
        }
    }

}