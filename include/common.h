//
// Created by Daniël de Vries on 26/11/2018.
//

#ifndef XROTOR_NOGRAPHICS_CPP_COMMON_H
#define XROTOR_NOGRAPHICS_CPP_COMMON_H

#include <string>
using namespace std;

namespace common {
    const int ix = 100;                 /**< max number of radial prop stations */
    const int ipx = ix + 1;
    const int nparx = 12;               /**< max number of stored cases */
    const int icasx = 100;              /**< number of case parameters stored */
    const int iwx = 200;                /**< dimension of work arrays */

    const int nax = 20;                 /**< max number of aerodynamic sections defined */
    const int ndx = 14;                 /**< number of aerodynamic parameters defined for each section */

    const int iq = ix + 5;
    const int jx = (iq * 3) / 2 + 1;

    const double pi = 3.14159265358979323846264338327950288;

    struct context {
        bool conv, greek, terse, vrtx, fast, free, duct, lstruc,
             dest, desp, stall, ldesini, loprint,
             lrotor, lvnorm, lpwrvar, legend;

        double mcrit;

        string savil, fname;
        string name;

        double rho, rmu, vso, vel, rad, gee, alt;

        int ii, iinf, incr, nn, nblds, ixspac,
            niterd, nitera;

        float version;
        double dt;

        int iaero[ix];

        double dbeta,
               xi0, xitip, xinf,
               xpitch, rake;
        double ch[ix], beta[ix], beta0[ix], t[ix],
               xi[ix], dxi[ix],
               xv[ix];

        int nadd;

        double ubody[ix];
        double urduct;

        int naero;
        double xiaero[nax];
        double aerodata[ndx][nax];

        int ncase, kcase, iwtyp;

        double adv, adw, adwfctr,
                rms, rlx, effinv,
                tspec, pspec, qspec,
                ttot, ptot, qtot,
                tinv, pinv, twak, pwak, tvis, pvis,
                gresmx, fresmx, aresmx;

        double w0[iwx], w1[iwx], w2[iwx], w3[iwx], w4[iwx],
               w5[iwx], w6[iwx], w7[iwx], w8[iwx], w9[iwx];

        int npwrvar;

        double xw0;

    };
}






#endif //XROTOR_NOGRAPHICS_CPP_COMMON_H
