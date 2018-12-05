//
// Created by daniel.devries on 11/30/2018.
//

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <iostream>

#include <spline.h>
#include <userio.h>
#include <xaero.h>
#include <xbend.h>
#include <xio.h>
#include <xnoise.h>
#include <xoper.h>
#include <xrotor.h>
using common::fprintf;

namespace xrotor {

    /**
     * Interactive Design and Analysis Program
     *     for Free-tip and Ducted Rotors
     */
    void XROTOR() {
        common::context ctxt;
        ctxt.VERSION = 7.55;

        cout << endl;
        cout << " =========================" << endl;
        cout << "    XROTOR Version " << ctxt.VERSION << endl;
        cout << " =========================" << endl;
        cout << endl;

        INIT(ctxt);

        string command, comarg;
        bool showHelp = true;
        while (true) {
            if (showHelp) {
                cout << "   QUIT   Exit program" << endl;
                cout << "  .OPER   Calculate off-design operating points" << endl;
                cout << "  .BEND   Calculate structural loads and deflections" << endl;
                cout << "  .NOIS   Calculate and plot acoustic signature" << endl;
                cout << "   LOAD f Read rotor from restart file" << endl;
                cout << "   DISP   Display current design point" << endl << endl;
                showHelp = false;
            }

            userio::ASKC(" XROTOR", command, comarg);

            ctxt.GREEK = true;
            if (command == "    ") continue;
            if (command == "?   ") {
                showHelp = true;
                continue;
            }
            if (command == "QUIT") return;
            if (command == "OPER") xoper::OPER(ctxt);
            if (command == "BEND") xbend::bend(ctxt);
            if (command == "LOAD") xio::LOAD(ctxt, comarg);
            if (command == "NOIS") xnoise::noise(ctxt);
            if (command == "DISP") OUTPUT(ctxt, cout);
            else if (ctxt.GREEK) {
                cout << command << " command not recognized.  Type a \"?\" for list" << endl;
            }
        }
    }

    /**
     * Initialize everything.
     *
     * @param ctxt
     */
    void INIT(common::context &ctxt) {
        ctxt.GREEK = false;

        // XROTOR defaults
        ctxt.URDUCT = 1.0;

        SETDEF(ctxt);

        if (ctxt.DUCT) {
            cout << "Aprop/Aexit initialized to 1.0" << endl;
            ctxt.URDUCT = 1.0;
        }

        ctxt.XINF = 3.0;                             // r/R at which BC at infinity is applied
        ctxt.NN = 32;                                // number of perturbation potential harmonics
        ctxt.IINF = ctxt.II + ctxt.II / 2;     // number of discrete potential harmonic stations
        ctxt.CONV = false;                           // operating point solution existence flag
        ctxt.LSTRUCT = false;                         // indicates if structural properties are available

        ctxt.NAME = " ";
        ctxt.SAVIL = " ";

        // acceleration due to gravity for scaling centrifugal blade tension (m/s^2)
        ctxt.GEE = 9.81;

        // ADW factor (multiplies TINV/PINV in ADW calculation)
        ctxt.ADWFCTR = 1.0;

        if (ctxt.II   > common::IX) throw runtime_error("Array overflow.  IX too small");
        if (ctxt.IINF > common::JX) throw runtime_error("Array overflow.  JX too small");

        // actual-rotor radius is always 1 (non-dimensionalized with itself)
        ctxt.XITIP = 1.0;

        // default nacelle, wake perturbation velocities (non-existent)
        for (double &i : ctxt.UBODY) i = 0;

        // no slipstream velocity profiles
        ctxt.NADD = 0;

        // number of defined cases
        ctxt.NCASE = 0;
        ctxt.KCASE = 0;

        // max number of iterations for design, analysis
        ctxt.NITERD = 40;
        ctxt.NITERA = 40;

        // do not initialize rotor at each design cycle
        ctxt.LDESINI = false;

        // do initialize rotor at each design cycle
        ctxt.LOPRINI = true;

        // no engine load line to start
        ctxt.LPWRVAR = false;
        ctxt.NPWRVAR = 0;
        
        // no rotor yet
        ctxt.LROTOR = false;
        for (int &i : ctxt.IAERO) i = 0;
    }

    /**
     * Hard-wired start-up defaults
     *
     * @param ctxt
     */
    void SETDEF(common::context &ctxt) {
        ctxt.RAKE = 0.0;

        ctxt.VEL = 1.0;
        ctxt.ALT = 0.0;
        ATMO(ctxt.ALT, ctxt.VSO, ctxt.RHO, ctxt.RMU); // sea level atmospheric conditions

        // install data into aero section #1
        xaero::PUTAERO(ctxt,
                       1, 0, 0,                    // NAERO, xisect, a0
                       1.5, -0.5, 6.28, 0.1, 0.1,  // clmax, clmin, dclda, dclda_stall, dcl_stall
                       0.013, 0.5, 0.004,          // cdmin, cldmin, cdccdl2
                       -0.1, 0.8,                  // cmcon, mcrit
                       200000, -0.4);              // reref, rexp
        for (int &i : ctxt.IAERO) i = 1;

        ctxt.XPITCH = 0.3;           // x/c location of pitch axis

        ctxt.II = 30;                // number of radial stations
        ctxt.INCR = 2;               // radial station increment for terminal OUTPUT
        ctxt.IXSPAC = 2;             // r/R spacing flag

        ctxt.VRTX = false;           // Vortex Wake (true)           / Graded Momentum (false) flag
        ctxt.FAST = false;           // Graded momentum (true)       / Potential Formulation (false) flag
        ctxt.FREE = true;            // Self-deforming wake (true)   / Rigid-wake (false) flag
        ctxt.DUCT = false;           // Ducted (true)                / Free-tip (false) flag

        ctxt.TERSE = false;          // TERSE-OUTPUT flag

        ctxt.LVNORM = true;          // flight speed used for normalization
    }

    /**
     * Calculate atmospheric properties at a given altitude.
     *
     * Returns speed of sound (VSO) in m/s, density (RHO)
     * in kg/m^3, and dynamic viscosity (RMU) in kg/m-s
     * of standard atmosphere at specified altitude ALSPEC
     * (in kilometers).  If ALSPEC=-1, water properties
     * at 15 Celsius are returned.
     *
     * Reference:  "U.S. Standard Atmosphere", NOAA.
     *
     * @param alspec        altitude in km
     * @param vsoalt        speed of sound in m/s
     * @param rhoalt        density in kg/m^3
     * @param rmualt        dynamic viscosity in kg/(m*s)
     */
    void ATMO(double alspec, double &vsoalt, double &rhoalt, double &rmualt) {
        const int n = 44;
        const bool first = true;
        const double alt[] = {  0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,
                               10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
                               20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0,
                               30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0,
                               40.0, 45.0, 60.0, 75.0};
        const double vso[] = {340.0,336.0,332.0,329.0,325.0,320.0,316.0,312.0,308.0,304.0,
                              299.0,295.0,295.0,295.0,295.0,295.0,295.0,295.0,295.0,295.0,
                              295.0,295.8,296.4,297.1,297.8,298.5,299.1,299.8,300.5,301.1,
                              301.8,302.5,303.1,305.0,306.8,308.7,310.5,312.3,314.0,316.0,
                              318.0,355.0,372.0,325.0};
        const double rho[] = {1.226,1.112,1.007,0.909,0.820,0.737,0.660,0.589,0.526,0.467,
                              0.413,0.364,0.311,0.265,0.227,0.194,0.163,0.141,0.121,0.103,
                              .0880,.0749,.0637,.0543,.0463,.0395,.0338,.0288,.0246,.0210,
                              .0180,.0154,.0132,.0113,.0096,.0082,.0070,.0060,.0052,.0044,
                              0.004,0.002,3.9E-4,8.0E-5};
        const double rmu[] = {1.780,1.749,1.717,1.684,1.652,1.619,1.586,1.552,1.517,1.482,
                              1.447,1.418,1.418,1.418,1.418,1.418,1.418,1.418,1.418,1.418,
                              1.418,1.427,1.433,1.438,1.444,1.449,1.454,1.460,1.465,1.471,
                              1.476,1.481,1.487,1.502,1.512,1.532,1.546,1.561,1.580,1.600,
                              1.700,1.912,2.047,1.667};

        // special case: Water at STP
        if (alspec == -1.0) {
            vsoalt = 1500;
            rhoalt = 1000;
            rmualt = 1.15e-3;
            cout << "                              o        " << endl;
            cout << "ATMO: You are underwater at 15  Celsius" << endl;
            return;
        }

        // linearly interpolate quantities from tabulated values
        for (int i = 1; i < n; i++) {
            if (alspec <= alt[i]) {
                double dalt = alt[i] - alt[i-1];
                double dvso = vso[i] - vso[i-1];
                double drho = rho[i] - rho[i-1];
                double drmu = rmu[i] - rmu[i-1];

                double alfrac = (alspec - alt[i-1]) / dalt;

                vsoalt = vso[i-1] + dvso * alfrac;
                rhoalt = rho[i-1] + drho * alfrac;
                rmualt = rmu[i-1] + drmu * alfrac;
                rmualt *= 1.0e-5;

                return;
            }
        }

        if (alspec > alt[n-1]) {
            cout << endl;
            cout << "ATMO: You''RE in low earth orbit.  Good luck." << endl;
            vsoalt = vso[n-1];
            rhoalt = vso[n-1];
            rmualt = rmu[n-1] * 1.0e-5;
        }
    }

    void FLOSHO(ostream &os, double vso, double rho, double rmu) {
        const double R = 287.0;
        const double gam = 1.4;
        double rnu = rmu / rho;
        double p = rho * pow(vso, 2) / gam;
        double t = p / (rho * R);
        cout << " Speed of sound (m/s): " << vso << endl;
        cout << " Density   (kg/m^3)  : " << rho << endl;
        cout << " Viscosity (kg/m-s)  : " << rmu << endl;
        cout << " Kin. Visc. (m^2/s)  : " << rnu << endl;
        cout << " Air pressure (Pa)   : " << p << endl;
        cout << " Air temperature (K) : " << t << endl;
    }

    /**
     * Re-initialize advance ratio and gammas.
     *
     * @param ctxt
     */
    void REINIT(common::context &ctxt) {
        // estimate reasonable advance ratio to start iterative routines
        int is = ctxt.II / 2 + 1;
        // HHY had to set A0 to 0.0 as A0 is now section property
        double a0 = 0;
        double ang = ctxt.BETA[is] - a0;

        double rpm = ctxt.VEL / (ctxt.RAD * ctxt.ADV * common::PI / 30.);

        double adv0 = ctxt.XI[is] * sin(ang) / cos(ang);
        double rpm0 = ctxt.VEL / (ctxt.RAD * adv0 * common::PI / 30.);

        userio::ASKR("Enter initialization RPM?", rpm);

        ctxt.ADV = ctxt.VEL / (rpm * ctxt.RAD * common::PI / 30.);
        ctxt.ADV = max(0.1, ctxt.ADV);
        ctxt.ADW = ctxt.ADV;

        // Set the blade angle back to reference angle
        bool yes;
        userio::ASKL("Restore blade angles to original?", yes);
        if (yes) copy(ctxt.BETA, ctxt.BETA + ctxt.II, ctxt.BETA0);

        // calculate current operating point
        xoper::APER(ctxt, 4, 2, true);
        if (ctxt.CONV) OUTPUT(ctxt, cout);
    }

    /**
     * Fill stretched radial coordinate array X (and XV).
     *
     * @param ctxt
     */
    void SETX(common::context &ctxt) {
        ctxt.DT = 0.5 * common::PI / float(ctxt.II);
        double xm = ctxt.XI0;
        ctxt.XV[0] = ctxt.XI0;

        double tp, xp;
        for (int i = 0; i < ctxt.II; i++) {
            ctxt.T[i] = ctxt.DT * (float(i) - 0.5);
            tp           = ctxt.DT * float(i);

            if (ctxt.IXSPAC == 2) {
                // Usual sine stretching, adjusted for nonzero root radius
                ctxt.XI[0] = sqrt(
                    ctxt.XITIP * pow(sin(ctxt.T[i]), 2) + pow(ctxt.XI0 * cos(ctxt.T[i]), 2));
                xp           = sqrt(
                    ctxt.XITIP * pow(sin(tp          ), 2) + pow(ctxt.XI0 * cos(tp          ), 2));
            } else {
                // Cosine stretching for more root resolution (also in TINVRT)
                ctxt.XI[0] = 0.5 * (1.0 - cos(2.0 * ctxt.T[i])) * (ctxt.XITIP - ctxt.XI0) + ctxt.XI0;
                xp            = 0.5 * (1.0 - cos(2.0 * tp          )) * (ctxt.XITIP - ctxt.XI0) + ctxt.XI0;
            }

            ctxt.XI[i] = (xp + xm) * 0.5;
            ctxt.DXI[i] = xp - xm;

            xm = xp;
            ctxt.XV[i+1] = xp;
        }
        ctxt.XV[ctxt.II + 1] = ctxt.XITIP;
    }

    /**
     * Open a file for writing.
     *
     * @param ofs       output file stream
     * @param fname     name of the file
     */
    void OPFILE(ofstream &ofs, string &fname) {
        // get filename if it hasn'T been already specified
        if (fname[0] == ' ') userio::ASKS("Enter output filename", fname);

        // try to open the file
        ofs.open(fname, ofstream::in);
        if (ofs.is_open()) {
            // file exists... ask how to proceed
            string command, comarg, valid;
            valid = "OoAaNn";
            userio::ASKC("File " + fname + " exists. Overwrite / Append / New file ?", command, comarg);

            if (valid.find(command[0]) == -1) {
                //ask again if reply is invalid
                userio::ASKC(" O / A / N ?", command, comarg);

                if (valid.find(command[0]) == -1) {
                    // Still bad reply. Give up asking and just return
                    cout << "No action taken" << endl;
                    return;
                }
            }

            // at this point, file is open and reply is valid
            valid = "Oo";
            if (valid.find(command[0]) != -1) {
                ofs.close();
                ofs.open(fname, ofstream::out | ofstream::trunc);
                return;
            }
            valid = "Aa";
            if (valid.find(command[0]) != -1) {
                ofs.close();
                ofs.open(fname, ofstream::out | ofstream::app);
                return;
            }

            // new file... get filename from command argument, or ask if not supplied
            fname = comarg;
            if (fname.empty() or fname.find_first_not_of(' ') == -1) userio::ASKS("Enter output filename", fname);
        } else {
            ofs.open(fname, ofstream::out | ofstream::trunc);
            if (!ofs.is_open()) cout << "Bad filename." << endl;
        }
    }

    /**
     * Dump everything to the FILE object identifying an output stream.
     *
     * @param ctxt
     * @param pFile     pointer to FILE object identifying output stream
     */
    void OUTPUT(common::context &ctxt, ostream &os) {
        int iadd = 1;
        if (&os == &cout) iadd = ctxt.INCR;

        fprintf(os, "%75s\n", string(75, '=').c_str());
        if (!ctxt.CONV) fprintf(os, "********** NOT CONVERGED **********\n");

        bool lheli = false;
        double fac1 = ctxt.RHO * pow(ctxt.VEL, 2) * pow(ctxt.RAD, 2);
        double fac2 = fac1 * ctxt.VEL;

        // dimensional thrust, power, torque, rpm
        double tdim = ctxt.TTOT * fac1;
        double qdim = ctxt.QTOT * fac1 * ctxt.RAD;
        double pdim = ctxt.PTOT * fac2 * ctxt.RAD;

        double tvdim = ctxt.TVIS * fac1;
        double pvdim = ctxt.PVIS * fac2;

        double efftot = ctxt.TTOT / ctxt.PTOT;
        double rpm = ctxt.VEL / (ctxt.RAD * ctxt.ADV * common::PI / 30.);
        double dia = 2.0 * ctxt.RAD;

        // Nacelle (or body) thrust is difference between thrust on
        // equivalent prop and real prop
        double tnacel = (ctxt.TWAK - ctxt.TINV) * fac1;

        // blade solidity
        spline::SPLINE(ctxt.CH, ctxt.W1, ctxt.XI, ctxt.II);
        double ch34 = spline::SEVAL(0.75, ctxt.CH, ctxt.W1, ctxt.XI, ctxt.II);
        double sigma = float(ctxt.NBLDS) * ch34 / common::PI;

        // standard coefficients based on forward speed
        double tc = tdim / (0.5 * fac1);
        double pc = pdim / (0.5 * fac2);

        // standard coefficients based on rotational speed
        double en = rpm / 60.;
        double ct = tdim / (ctxt.RHO * pow(en, 2) * pow(dia, 4));
        double cp = pdim / (ctxt.RHO * pow(en, 3) * pow(dia, 5));

        // induced efficiency (including nacelle thrust effect)
        double effind = ctxt.TWAK / ctxt.PWAK;

        // ideal (actuator disk) efficiency
        double tclim = max(-1.0, tc);
        double eideal = 2.0 / (1.0 + sqrt(tclim + 1.0));

        double cth, cph, ctos, fom;
        // define low advance ratio (helicopter?) related data
        if (ctxt.ADV < 0.1) {
            spline::SPLINE(ctxt.CH, ctxt.W1, ctxt.XI, ctxt.II);
            cth = ct / 7.7516;
            cph = cp / 24.352;
            ctos = cth / sigma;
            fom = 0.7979 * pow(abs(ct), 1.5) / cp;
            lheli = true;
        }

        if (ctxt.DUCT) {
            switch (ctxt.IWTYP) {
                case 1:
                case 3: fprintf(os, " Ducted Graded Mom. Formulation Solution:  ");
                case 2: fprintf(os, " Ducted Potential Formulation Solution:  ");
                default:;
            }
        } else {
            switch (ctxt.IWTYP) {
                case 1: fprintf(os, " Free Tip Graded Mom. Formulation Solution:  ");
                case 2: fprintf(os, " Free Tip Potential Formulation Solution:  ");
                case 3: fprintf(os, " Free Tip Vortex Wake Formulation Solution:  ");
                default:;
            }
        }
        fprintf(os, "%32s\n", ctxt.NAME.c_str());

        if (ctxt.NADD > 1) {
            fprintf(os, " (External slipstream present)%19sWake ADV. ratio:%11.5f\n", "", ctxt.ADW);
        } else if (ctxt.DUCT) {
            fprintf(os, " Vdisk/Vslip:%11.5f%25sWake ADV. ratio:%11.5f\n", ctxt.URDUCT, "", ctxt.ADW);
        } else {
            fprintf(os, "%50sWake ADV. ratio:%11.5f\n", "", ctxt.ADW);
        }

        if (ctxt.ADW < 0.5*ctxt.ADV)
            fprintf(os, " Reverse far-slipstream velocity implied. Interpret results carefully !");

        fprintf(os, " no. blades :%3i   %9sradius(m)  :%9.4F %4sadv. ratio: %11.5F\n", ctxt.NBLDS, "", ctxt.RAD, "", ctxt.ADV);
        fprintf(os, " thrust(N)  :%11.3G%4spower(W)   :%11.3G%3storque(N-m):%11.3G\n", tdim, "", pdim, "", qdim);
        fprintf(os, " Efficiency :%8.3F %6sspeed(m/s) :%9.3F %4srpm        :%11.3F\n", efftot, "", ctxt.VEL, "", rpm);
        fprintf(os, " Eff induced:%8.4F %6sEff ideal  :%9.4F %4sTcoef      :%11.4F\n", effind, "", eideal, "", tc);
        fprintf(os, " Tnacel(N)  :%11.4F%4shub RAD.(m):%9.4F %4sdisp. RAD. :%10.4F\n", tnacel, "", ctxt.XI0 * ctxt.RAD, "", ctxt.XW0 * ctxt.RAD);
        fprintf(os, " Tvisc(N)   :%11.4F%4sPvisc(W)   :%11.3G\n",                      tvdim, "", pvdim);
        fprintf(os, " RHO(kg/m3) :%10.5F%5sVsound(m/s):%9.3F %4smu(kg/m-s) :%11.4E\n", ctxt.RHO, "", ctxt.VSO, "", ctxt.RMU);
        fprintf(os, " %75s\n", string(75, '-').c_str());

        // low advance ratio (helicopter?) data
        if (lheli) {
            fprintf(os, "Helicopter: Sigma:%11.5F  CTh/s:%11.5F  FOM:%11.5F", sigma, ctos, fom);
        } else {
            fprintf(os, " Sigma:%11.5F", sigma);
        }

        // coefficients based on rotational speed
        fprintf(os, "%12s    Ct:%11.5F     Cp:%11.5F    J:%11.5F", "", ct, cp, ctxt.ADV * common::PI);
        // coefficients based on forward speed
        fprintf(os, "%12s    Tc:%11.5F     Pc:%11.5F  ADV:%11.5F", "", tc, pc, ctxt.ADV);

        if (ctxt.TERSE) return;

        // find maximum RE on blade
        double remax = 0.0;
        for (int i = 0; i < ctxt.II; i++) {
            remax = max(ctxt.RE[i], remax);
        }
        double reexp = 1.0;
        if (remax >= 1.0e6) {
            reexp = 6.0;
        } else if (remax >= 1.0e3) {
            reexp = 3.0;
        }

        if (reexp == 1.0) {
            fprintf(os, "\n  i  r/R    c/R  BETA(deg)   CL      Cd    RE    Mach   effi  EFFP  na.u/U");
        } else {
            fprintf(os, "\n  i  r/R   c/R  BETA(deg)  CL     Cd    REx10^%1i Mach   effi  EFFP  na.u/U", (int)reexp);
        }

        double wa, wt;
        double vw, vaw, utotw, cw, sw, effi;
        double utot, vt, vt_adw, va, va_adw, vd, vd_adw, ci, ci_adv, ci_vt,
               si, si_va, w, w_adv, w_vt, w_va, phi, p_adv, p_vt, p_va;
        double mach, bdeg, xre;
        char schar[1] = {' '};
        for (int i = 0; i < ctxt.II; i += iadd) {
            // use equivalent prop to define local efficiency
            UVADD(ctxt, ctxt.XI[i], wa, wt);
            vw = ctxt.VWAK[i];
            vaw = vw * ctxt.XW[i] / ctxt.ADW;
            // Freestream velocity component on equiv prop
            utotw = ctxt.URDUCT;
            cw = ctxt.XI[i] / ctxt.ADV - wt - vw;
            sw = utotw                       + wa - vaw;
            effi = (cw / sw)                 * ctxt.ADV / ctxt.XW[i];

            // use real prop to define Mach number
            xoper::CSCALC(ctxt,
                          i, utot, wa, wt,
                          vt, vt_adw,
                          va, va_adw,
                          vd, vd_adw,
                          ci, ci_adv, ci_vt,
                          si, si_va,
                          w, w_adv, w_vt, w_va,
                          phi, p_adv, p_vt, p_va);

            mach = w * ctxt.VEL / ctxt.VSO;

            bdeg = ctxt.BETA[i] * 180. / common::PI;
            xre = ctxt.RE[i] / (pow(10., reexp));

            schar[0] = ' ';
            if (ctxt.STALL[i]) schar[0] = 'S';

            fprintf(os, " %2i%6.3F%7.4F%7.2F%7.3F %1s%7.4F %6.2F %6.3F %6.3F%6.3F%8.3F\n",
                i, ctxt.XI[i], ctxt.CH[i], bdeg, ctxt.CL[i], schar, ctxt.CD[i], xre, mach,
                effi, ctxt.EFFP[i], ctxt.UBODY[i]);
        }
    }

    /**
     * @param ctxt
     * @param xiw
     * @param wa
     * @param wt
     */
    void UVADD(common::context &ctxt, double xiw, double &wa, double &wt) {
        wa = 0;
        wt = 0;

        if (ctxt.NADD <= 1) return;

        double rdim = xiw * ctxt.RAD;
        if (rdim >= ctxt.RADD[0] and rdim <= ctxt.RADD[ctxt.NADD-1]) {
            wa = spline::SEVAL(rdim, ctxt.UADD, ctxt.UADDR, ctxt.RADD, ctxt.NADD) / ctxt.VEL;
            wt = spline::SEVAL(rdim, ctxt.VADD, ctxt.VADDR, ctxt.RADD, ctxt.NADD) / ctxt.VEL;
        }
    }
}