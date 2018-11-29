//
// Created by Daniël de Vries on 2018-11-28.
//
#include <cmath>        // std::sin, std::cos, std::tan, std::sqrt
#include <iostream>     // std::cout, std::endl
#include "common.h"
#include "vortex.h"
using namespace std;

namespace vortex {

    /**
     * Calculate "Vortex Momentum" Gamma-swirl influence coefficients.
     *
     * @param ii        number radial points on blade (circulation stations)
     * @param nblds     number of blades
     * @param lduct     true for duct outer BC
     * @param rake      blade rake angle from y-axis in xy-plane
     * @param xi        r/R control point coordinate array
     * @param xv        r/R vortex leg    coordinate array
     * @param gam       circulation array
     * @param adw       wake advance ratio V/wR
     * @param vind_gam  sensitivity of swirl velocity to circulation (calculated)
     *                  @note must of of size (3, ii, ii)
     * @param vind_adw  sensitivity of swirl velocity to wake advance ratio (calculated)
     *                  @note must be of size (3, ii)
     */
    void vrtxc0(unsigned long ii, unsigned nblds, bool lduct, double rake,
                const vec &xi, const vec &xv, const vec &gam, double adw,
                Cube &vind_gam, Matrix &vind_adw) {
        auto blds = (double)nblds;
        double pi = common::pi;

        double xi0 = xv[0];
        double xitip = xv[ii + 1];
        double tanrak = tan(rake);

        vind_gam.zeros();
        vind_adw.zeros();

        // Set up variable theta spacing for near, intermediate and far field
        double dth1  = pi / 60.;
        double rad1  = 2.0;
        double thet1 = rad1 / adw;
        double dth2  = pi / 20.;
        double rad2  = 4.0;
        double thet2 = rad2 / adw;
        double dth3  = pi / 8.;
        double rad3  = 50.0;
        double thet3 = rad3 / adw;

        double ddth1 = (dth2 - dth1) / (2.0 * (thet1        ) / (dth2 + dth3) - 1.0);
        double ddth2 = (dth3 - dth2) / (2.0 * (thet2 - thet1) / (dth3 + dth2) - 1.0);

        const int ntdim = 5000;
        vec thetspc(ntdim);

        unsigned long nthet;
        double thet = 0.0;
        double dth = dth1;
        for (unsigned long i = 0; true; i++) {
            if (i < ntdim) {
                thetspc[i] = thet;
            } else {
                // TODO: This approach deviates from the FORTRAN code. Check that this actually works.
                thetspc.push_back(thet);
            }
            thet += dth;

            if (thet < thet1) {
                dth += ddth1;
            } else if (thet < thet2) {
                dth += ddth2;
            } else if (thet < thet3) {
                dth = dth3;
            } else {
                nthet = i-1;
                break;
            }
        }
        if (nthet > ntdim) {
            cout << "The number of vortex segments is very large" << endl;
        }
        // TODO: In the FORTRAN code, the number of segments is capped at ntdim. Check that this works like this too.

        if (lduct) {
            // use simple mean swirl to get swirl at blade
            for (unsigned long i = 0; i < ii; i++) {
                vind_gam(2, i, i) =  blds / (4.0 * pi * xi[i]);
                vind_gam(0, i, i) =  vind_gam(2, i, i) * xi[i]  / adw;
                vind_adw(0, i)    = -vind_gam(0, i, i) * gam[i] / adw;
            }
        } else {
            // Do a discrete vortex integration of slipstream vortices
            double dtbld = 2.0 * pi / blds;
            cout << "Vortex points/radial station = " << nblds * nthet;

            double r0x, r0y, r0z,
                   vsum[3], vadw[3], rv, xxv, thetoff,
                   r1x, r1y, r1z, r1_adw,
                   r2x, r2y, r2z, r2_adw, a[3], b[3], uvw[3];
            Matrix uvw_a(3, 3), uvw_b(3, 3);
            // velocity influences for point - R0
            for (unsigned long i = 0; i < ii; i++) {
                r0x = xi[i] * tanrak;
                r0y = xi[i];
                r0z = 0.0;

                // For each vortex trailing leg (II+1 legs starting at XV(J))
                for (unsigned long j = 0; j < ii+1; j++) {
                    vsum[0] = 0; vsum[1] = 0; vsum[2] = 0;
                    vadw[0] = 0; vadw[1] = 0; vadw[2] = 0;

                    rv  = xv[j];
                    xxv = rv * tanrak;

                    // For each blade
                    thetoff = 0.0;
                    for (unsigned n = 0; n < nblds; n++) {
                        // For angles around helix to the far-field
                        thet1 = thetspc[0];
                        r1x = xxv + thet1 * xitip * adw;
                        r1y = rv * cos(thet1 + thetoff);
                        r1z = rv * sin(thet1 + thetoff);
                        r1_adw = thet1 * xitip;

                        for (unsigned long l = 0; l < nthet - 1; l++) {
                            thet2 = thetspc[l + 1];
                            r2x = xxv + thet2 * xitip * adw;
                            r2y = rv * cos(thet2 + thetoff);
                            r2z = rv * sin(thet2 + thetoff);
                            r2_adw = thet2 * xitip;

                            a[0] = r1x - r0x;
                            a[1] = r1y - r0y;
                            a[2] = r1z - r0z;
                            b[0] = r2x - r0x;
                            b[1] = r2y - r0y;
                            b[2] = r2z - r0z;
                            vortex::vorsegvel(a, b, uvw, uvw_a, uvw_b);

                            vsum[0] += uvw[0];
                            vsum[1] += uvw[1];
                            vsum[2] += uvw[2];

                            vadw[0] += uvw_a(0, 0) * r1_adw
                                    +  uvw_b(0, 0) * r2_adw;
                            vadw[1] += uvw_a(1, 1) * r1_adw
                                    +  uvw_b(1, 1) * r2_adw;
                            vadw[2] += uvw_a(2, 2) * r1_adw
                                    +  uvw_b(2, 2) * r2_adw;

                            thet1 = thet2;
                            r1x = r2x;
                            r1y = r2y;
                            r1z = r2z;
                            r1_adw = r2_adw;
                        }

                        thetoff += dtbld;
                    }

                    vsum[2] = -vsum[2];
                    vadw[2] = -vadw[2];
                    // alternate + and - influence for each vortex line into velocity influence matrix

                    // Open wake, interdigitate all vortex lines
                    if (j < ii) {
                        vind_gam(0, i, j) = -vsum[0];
                        vind_gam(1, i, j) = -vsum[1];
                        vind_gam(2, i, j) = -vsum[2];
                        vind_adw(0, i) -= vadw[0] * gam[j];
                        vind_adw(1, i) -= vadw[1] * gam[j];
                        vind_adw(2, i) -= vadw[2] * gam[j];
                    }
                    if (j > 0) {
                        vind_gam(0, i, j-1) -= vsum[0];
                        vind_gam(1, i, j-1) -= vsum[1];
                        vind_gam(2, i, j-1) -= vsum[2];
                        vind_adw(0, i) += vadw[0] * gam[j-1];
                        vind_adw(1, i) += vadw[1] * gam[j-1];
                        vind_adw(2, i) += vadw[2] * gam[j-1];
                    }
                }
            }
        }

    }

    /**
     * Calculate the velocity induced by a vortex segment of unit strength, with no core radius.
     *
     * The point where the velocity is calculated is at (0, 0, 0).
     *
     * Positive circulation is by righthand rule from A to B.
     *
     * @param a     coordinates of vertex #1 of the vortex
     * @param b     coordinates of vertex #2 of the vortex
     * @param uvw   incuded velocity (calculated)
     *              @note must be of size (3)
     * @param uvw_a dUVW/dA sensitivity (calculated)
     *              @note must be of size (3, 3)
     * @param uvw_b dUVW/dB sensitivity (calculated)
     *              @note must be of size (3, 3)
     */
    void vorsegvel(const double a[3], const double b[3], double uvw[3], Matrix &uvw_a, Matrix &uvw_b) {
        double asq = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
        double bsq = b[0]*b[0] + b[1]*b[1] + b[2]*b[2];

        double amag = sqrt(asq);
        double bmag = sqrt(bsq);

        uvw[0] = 0; uvw[1] = 0; uvw[2] = 0;
        uvw_a.zeros();
        uvw_b.zeros();

        // contribution from the vortex leg
        if (amag * bmag != 0) {
            double axb[3] = {a[1]*b[2] - a[2]*b[2],
                             a[2]*b[0] - a[0]*b[2],
                             a[0]*b[1] - a[1]*b[0]};

            double axb_a[3][3] = {
                    { 0.0 ,  b[2], -b[1]},
                    {-b[2],  0.0 ,  b[0]},
                    { b[1], -b[0],  0.0 }
            };

            double axb_b[3][3] = {
                    { 0.0 ,  a[2], -a[1]},
                    {-a[2],  0.0 ,  a[0]},
                    { a[1], -a[0],  0.0 }
            };

            double adb = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

            double den     =       amag * bmag + adb;
            double den_asq = 0.5 * bmag / amag;
            double den_bsq = 0.5 * amag / bmag;

            double t = (1.0 / amag + 1.0 / bmag) / den;

            double t_adb = -t / den;
            double t_asq = -t / den * den_asq - 0.5 / (den * amag * asq);
            double t_bsq = -t / den * den_bsq - 0.5 / (den * bmag * bsq);

            for (unsigned k = 0; k < 2; k++) {
                uvw[k] = axb[k] * t;

                for (unsigned l = 0; l < 2; l++) {
                    uvw_a(k, l) = (axb[k] * t_asq) * (a[l] * 2.0)
                                + (axb[k] * t_adb) * b[l]
                                + (axb_a[k][l] * t);
                    uvw_b(k, l) = (axb[k] * t_bsq) * (b[l] * 2.0)
                                  + (axb[k] * t_adb) * a[l]
                                  + (axb_b[k][l] * t);
                }
            }
        }

        double pi4 = 4.0 * common::pi;
        for (unsigned k = 0; k < 2; k++) {
            uvw[k] /= pi4;

            for (unsigned l = 0; l < 2; l++) {
                uvw_a(k, l) /= pi4;
                uvw_b(k, l) /= pi4;
            }
        }
    }

}
