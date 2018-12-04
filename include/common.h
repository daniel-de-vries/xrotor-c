//
// Created by Daniël de Vries on 26/11/2018.
//

#ifndef XROTOR_NOGRAPHICS_CPP_COMMON_H
#define XROTOR_NOGRAPHICS_CPP_COMMON_H

#include <fstream>
#include <cstdio>
#include <cstdarg>
#include <string>
using namespace std;

namespace common {
    const int IX = 100;                 /**< max number of radial prop stations */
    const int IPX = IX + 1;
    const int NPARX = 12;               /**< max number of stored cases */
    const int ICASX = 100;              /**< number of case parameters stored */
    const int IWX = 200;                /**< dimension of work arrays */

    const int NAX = 20;                 /**< max number of aerodynamic sections defined */
    const int NDX = 14;                 /**< number of aerodynamic parameters defined for each section */

    const int JQ = IX + 5;
    const int JX = (JQ * 3) / 2 + 1;

    const double PI = 3.14159265358979323846264338327950288;

    struct context {
        double MCRIT;

        bool CONV, GREEK, TERSE, VRTX, FAST, FREE, DUCT, LSTRUCT,
            DEST, DESP, LDESINI, LOPRINI,
            LROTOR, LVNORM, LPWRVAR, STALL[IX], LEGEND;

        string SAVIL, FNAME, NAME;

        double RHO, RMU, VSO, VEL, RAD, GEE, ALT;

        int II, IINF, INCR, NN, NBLDS, IXSPAC,
            NITERD, NITERA;

        float VERSION;
        double DT;

        int IAERO[IX];

        double CH[IX], BETA[IX], BETA0[IX], T[IX], DBETA,
               XI[IX], DXI[IX], XI0, XITIP, XINF,
               XPITCH, XV[IX], RAKE;

        int NADD;
        double RADD[IX],
               UADD[IX], VADD[IX],
               UADDR[IX], VADDR[IX],
               UBODY[IX], VBODY[IX], URDUCT;

        double CL[IX], CD[IX], CM[IX],
               RE[IX], EFFP[IX], GAM[IX],
               DTII[IX], DPII[IX],
               DTVI[IX], DPVI[IX],
               DTWI[IX], DPWI[IX];

        int NAERO;
        double XIAERO[NAX], AERODATA[NDX][NAX];

        int NCASE, KCASE, IWTYP;

        double ADV, ADW, ADWFCTR,
               RMS, RLX, EFFINV,
               TSPEC, PSPEC, QSPEC,
               TTOT, PTOT, QTOT,
               TINV, PINV, TWAK, PWAK, TVIS, PVIS,
               GRESMX, FRESMX, ARESMX;

        double W0[IWX], W1[IWX], W2[IWX], W3[IWX], W4[IWX],
               W5[IWX], W6[IWX], W7[IWX], W8[IWX], W9[IWX];

        int NPWRVAR;

        double VWAK[IX], VW_GAM[IX][IX], VW_ADW[IX], VW_ADV[IX],
               VIND[3][IX], VIND_GAM[3][IX][IX], VIND_ADW[3][IX];

        double XW0, XWTIP,
               XW[IX],   xw_gam[IX][IX],  xw_adw[IX][IX],  xw_adv[IX][IX],
               dwx[IX], dxw_gam[IX][IX], dxw_adw[IX][IX], dxw_adv[IX][IX];

    };

    /**
     * Write a formatted string to an output stream.
     *
     * @param os        output stream
     * @param format    format string
     * @param ...       arguments to the format string
     */
    void fprintf(ostream& os, const char* format, ...) {
        va_list argv;
        va_start(argv, format);
        char buff[1024];
        vsprintf(buff, format, argv);
        va_end(argv);
        os << buff;
    }
}






#endif //XROTOR_NOGRAPHICS_CPP_COMMON_H
