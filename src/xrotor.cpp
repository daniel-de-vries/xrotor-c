//
// Created by daniel.devries on 11/30/2018.
//

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <iostream>
#include <stdio.h>
#include "spline.h"
#include "userio.h"
#include "xaero.h"
#include "xbend.h"
#include "xio.h"
#include "xnoise.h"
#include "xoper.h"
#include "xrotor.h"


namespace xrotor {

    /**
     * Interactive Design and Analysis Program
     *     for Free-tip and Ducted Rotors
     */
    void xrotor() {
        common::context context;
        context.version = 7.55;

        cout << endl;
        cout << " =========================" << endl;
        cout << "    XROTOR Version " << context.version << endl;
        cout << " =========================" << endl;
        cout << endl;

        init(context);

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

            userio::askc(" XROTOR", command, comarg);

            context.greek = true;
            if (command == "    ") continue;
            if (command == "?   ") {
                showHelp = true;
                continue;
            }
            if (command == "QUIT") return;
            if (command == "OPER") xoper::oper(context);
            if (command == "BEND") xbend::bend(context);
            if (command == "LOAD") xio::load(context, comarg);
            if (command == "NOIS") xnoise::noise(context);
            if (command == "DISP") output(context, cout);
            else if (context.greek) {
                cout << command << " command not recognized.  Type a \"?\" for list" << endl;
            }
        }
    }

    /**
     * Initialize everything.
     *
     * @param context
     */
    void init(common::context& context) {
        context.greek = false;

        // XROTOR defaults
        context.urduct = 1.0;

        setdef(context);

        if (context.duct) {
            cout << "Aprop/Aexit initialized to 1.0" << endl;
            context.urduct = 1.0;
        }

        context.xinf = 3.0;                             // r/R at which BC at infinity is applied
        context.nn = 32;                                // number of perturbation potential harmonics
        context.iinf = context.ii + context.ii / 2;     // number of discrete potential harmonic stations
        context.conv = false;                           // operating point solution existence flag
        context.lstruc = false;                         // indicates if structural properties are available

        context.name = " ";
        context.savil = " ";

        // acceleration due to gravity for scaling centrifugal blade tension (m/s^2)
        context.gee = 9.81;

        // ADW factor (multiplies TINV/PINV in ADW calculation)
        context.adwfctr = 1.0;

        if (context.ii   > common::ix) throw runtime_error("Array overflow.  IX too small");
        if (context.iinf > common::jx) throw runtime_error("Array overflow.  JX too small");

        // actual-rotor radius is always 1 (non-dimensionalized with itself)
        context.xitip = 1.0;

        // default nacelle, wake perturbation velocities (non-existent)
        for (double &i : context.ubody) i = 0;

        // no slipstream velocity profiles
        context.nadd = 0;

        // number of defined cases
        context.ncase = 0;
        context.kcase = 0;

        // max number of iterations for design, analysis
        context.niterd = 40;
        context.nitera = 40;

        // do not initialize rotor at each design cycle
        context.ldesini = false;

        // do initialize rotor at each design cycle
        context.loprint = true;

        // no engine load line to start
        context.lpwrvar = false;
        context.npwrvar = 0;
        
        // no rotor yet
        context.lrotor = false;
        for (int &i : context.iaero) i = 0;
    }

    /**
     * Hard-wired start-up defaults
     *
     * @param context
     */
    void setdef(common::context& context) {
        context.rake = 0.0;

        context.vel = 1.0;
        context.alt = 0.0;
        atmo(context.alt, context.vso, context.rho, context.rmu); // sea level atmospheric conditions

        // install data into aero section #1
        xaero::putaero(context,
            1, 0, 0,                    // naero, xisect, a0
            1.5, -0.5, 6.28, 0.1, 0.1,  // clmax, clmin, dclda, dclda_stall, dcl_stall
            0.013, 0.5, 0.004,          // cdmin, cldmin, cdccdl2
            -0.1, 0.8,                  // cmcon, mcrit
            200000, -0.4);              // reref, rexp
        for (int &i : context.iaero) i = 1;

        context.xpitch = 0.3;           // x/c location of pitch axis

        context.ii = 30;                // number of radial stations
        context.incr = 2;               // radial station increment for terminal output
        context.ixspac = 2;             // r/R spacing flag

        context.vrtx = false;           // Vortex Wake (true)           / Graded Momentum (false) flag
        context.fast = false;           // Graded momentum (true)       / Potential Formulation (false) flag
        context.free = true;            // Self-deforming wake (true)   / Rigid-wake (false) flag
        context.duct = false;           // Ducted (true)                / Free-tip (false) flag

        context.terse = false;          // terse-output flag

        context.lvnorm = true;          // flight speed used for normalization
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
    void atmo(double alspec, double &vsoalt, double &rhoalt, double &rmualt) {
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
            cout << "ATMO: You''re in low earth orbit.  Good luck." << endl;
            vsoalt = vso[n-1];
            rhoalt = vso[n-1];
            rmualt = rmu[n-1] * 1.0e-5;
        }
    }

    void flosho(ostream &os, double vso, double rho, double rmu) {
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
     * @param context
     */
    void reinit(common::context &context) {
        // estimate reasonable advance ratio to start iterative routines
        int is = context.ii / 2 + 1;
        // HHY had to set A0 to 0.0 as A0 is now section property
        double a0 = 0;
        double ang = context.beta[is] - a0;

        double rpm = context.vel / (context.rad * context.adv * common::pi / 30.);

        double adv0 = context.xi[is] * sin(ang) / cos(ang);
        double rpm0 = context.vel / (context.rad * adv0 * common::pi / 30.);

        userio::askr("Enter initialization RPM?", rpm);

        context.adv = context.vel / (rpm * context.rad * common::pi / 30.);
        context.adv = max(0.1, context.adv);
        context.adw = context.adv;

        // Set the blade angle back to reference angle
        bool yes;
        userio::askl("Restore blade angles to original?", yes);
        if (yes) copy(context.beta, context.beta + context.ii, context.beta0);

        // calculate current operating point
        xoper::aper(context, 4, 2, true);
        if (context.conv) output(context, cout);
    }

    /**
     * Fill stretched radial coordinate array X (and XV).
     *
     * @param context
     */
    void setx(common::context &context) {
        context.dt = 0.5 * common::pi / float(context.ii);
        double xm = context.xi0;
        context.xv[0] = context.xi0;

        double tp, xp;
        for (int i = 0; i < context.ii; i++) {
            context.t[i] = context.dt * (float(i) - 0.5);
            tp           = context.dt * float(i);

            if (context.ixspac == 2) {
                // Usual sine stretching, adjusted for nonzero root radius
                context.xi[0] = sqrt(
                    context.xitip * pow(sin(context.t[i]), 2) + pow(context.xi0 * cos(context.t[i]), 2));
                xp           = sqrt(
                    context.xitip * pow(sin(tp          ), 2) + pow(context.xi0 * cos(tp          ), 2));
            } else {
                // Cosine stretching for more root resolution (also in TINVRT)
                context.xi[0] = 0.5 * (1.0 - cos(2.0 * context.t[i])) * (context.xitip - context.xi0) + context.xi0;
                xp            = 0.5 * (1.0 - cos(2.0 * tp          )) * (context.xitip - context.xi0) + context.xi0;
            }

            context.xi[i] = (xp + xm) * 0.5;
            context.dxi[i] = xp - xm;

            xm = xp;
            context.xv[i+1] = xp;
        }
        context.xv[context.ii + 1] = context.xitip;
    }

    /**
     * Open a file for writing.
     *
     * @param ofs       output file stream
     * @param fname     name of the file
     */
    void opfile(ofstream &ofs, string &fname) {
        // get filename if it hasn't been already specified
        if (fname[0] == ' ') userio::asks("Enter output filename", fname);

        // try to open the file
        ofs.open(fname, ofstream::in);
        if (ofs.is_open()) {
            // file exists... ask how to proceed
            string command, comarg, valid;
            valid = "OoAaNn";
            userio::askc("File " + fname + " exists. Overwrite / Append / New file ?", command, comarg);

            if (valid.find(command[0]) == -1) {
                //ask again if reply is invalid
                userio::askc(" O / A / N ?", command, comarg);

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
            if (fname.empty() or fname.find_first_not_of(' ') == -1) userio::asks("Enter output filename", fname);
        } else {
            ofs.open(fname, ofstream::out | ofstream::trunc);
            if (!ofs.is_open()) cout << "Bad filename." << endl;
        }
    }

    /**
     * Dump everything to the FILE object identifying an output stream.
     *
     * @param context
     * @param pFile     pointer to FILE object identifying output stream
     */
    void output(common::context &context, FILE* pFile) {
        int iadd = 1;
        if (pFile == stdout) iadd = context.incr;

        fprintf(pFile, "%75s\n", string(75, '=').c_str());
        if (!context.conv) fprintf(pFile, "********** NOT CONVERGED **********\n");

        bool lheli = false;
        double fac1 = context.rho * pow(context.vel, 2) * pow(context.rad, 2);
        double fac2 = fac1 * context.vel;

        // dimensional thrust, power, torque, rpm
        double tdim = context.ttot * fac1;
        double qdim = context.qtot * fac1 * context.rad;
        double pdim = context.ptot * fac2 * context.rad;

        double tvdim = context.tvis * fac1;
        double pvdim = context.pvis * fac2;

        double efftot = context.ttot / context.ptot;
        double rpm = context.vel / (context.rad * context.adv * common::pi / 30.);
        double dia = 2.0 * context.rad;

        // Nacelle (or body) thrust is difference between thrust on
        // equivalent prop and real prop
        double tnacel = (context.twak - context.tinv) * fac1;

        // blade solidity
        spline::spline(context.ch, context.w1, context.xi, context.ii);
        double ch34 = spline::seval(0.75, context.ch, context.w1, context.xi, context.ii);
        double sigma = float(context.nblds) * ch34 / common::pi;

        // standard coefficients based on forward speed
        double tc = tdim / (0.5 * fac1);
        double pc = pdim / (0.5 * fac2);

        // standard coefficients based on rotational speed
        double en = rpm / 60.;
        double ct = tdim / (context.rho * pow(en, 2) * pow(dia, 4));
        double cp = pdim / (context.rho * pow(en, 3) * pow(dia, 5));

        // induced efficiency (including nacelle thrust effect)
        double effind = context.twak / context.pwak;

        // ideal (actuator disk) efficiency
        double tclim = max(-1.0, tc);
        double eideal = 2.0 / (1.0 + sqrt(tclim + 1.0));

        double cth, cph, ctos, fom;
        // define low advance ratio (helicopter?) related data
        if (context.adv < 0.1) {
            spline::spline(context.ch, context.w1, context.xi, context.ii);
            cth = ct / 7.7516;
            cph = cp / 24.352;
            ctos = cth / sigma;
            fom = 0.7979 * pow(abs(ct), 1.5) / cp;
            lheli = true;
        }

        if (context.duct) {
            switch (context.iwtyp) {
                case 1:
                case 3: fprintf(pFile, " Ducted Graded Mom. Formulation Solution:  ");
                case 2: fprintf(pFile, " Ducted Potential Formulation Solution:  ");
                default:;
            }
        } else {
            switch (context.iwtyp) {
                case 1: fprintf(pFile, " Free Tip Graded Mom. Formulation Solution:  ");
                case 2: fprintf(pFile, " Free Tip Potential Formulation Solution:  ");
                case 3: fprintf(pFile, " Free Tip Vortex Wake Formulation Solution:  ");
                default:;
            }
        }
        fprintf(pFile, "%32s\n", context.name.c_str());

        if (context.nadd > 1) {
            fprintf(pFile, " (External slipstream present)%19sWake adv. ratio:%11.5f\n", "", context.adw);
        } else if (context.duct) {
            fprintf(pFile, " Vdisk/Vslip:%11.5f%25sWake adv. ratio:%11.5f\n", context.urduct, "", context.adw);
        } else {
            fprintf(pFile, "%50sWake adv. ratio:%11.5f\n", "", context.adw);
        }

        if (context.adw < 0.5*context.adv)
            fprintf(pFile, " Reverse far-slipstream velocity implied. Interpret results carefully !");

        fprintf(pFile, " no. blades :%3i   %9sradius(m)  :%9.4F %4sadv. ratio: %11.5F\n", context.nblds, "", context.rad, "", context.adv);
        fprintf(pFile, " thrust(N)  :%11.3G%4spower(W)   :%11.3G%3storque(N-m):%11.3G\n", tdim, "", pdim, "", qdim);
        fprintf(pFile, " Efficiency :%8.3F %6sspeed(m/s) :%9.3F %4srpm        :%11.3F\n", efftot, "", context.vel, "", rpm);
        fprintf(pFile, " Eff induced:%8.4F %6sEff ideal  :%9.4F %4sTcoef      :%11.4F\n", effind, "", eideal, "", tc);
        fprintf(pFile, " Tnacel(N)  :%11.4F%4shub rad.(m):%9.4F %4sdisp. rad. :%10.4F\n", tnacel, "", context.xi0 * context.rad, "", context.xw0 * context.rad);
        fprintf(pFile, " Tvisc(N)   :%11.4F%4sPvisc(W)   :%11.3G\n",                      tvdim, "", pvdim);
        fprintf(pFile, " rho(kg/m3) :%10.5F%5sVsound(m/s):%9.3F %4smu(kg/m-s) :%11.4E\n", context.rho, "", context.vso, "", context.rmu);
        fprintf(pFile, " %75s\n", string(75, '-').c_str());

        // low advance ratio (helicopter?) data
        if (lheli) {
            fprintf(pFile, "Helicopter: Sigma:%11.5F  CTh/s:%11.5F  FOM:%11.5F", sigma, ctos, fom);
        } else {
            fprintf(pFile, " Sigma:%11.5F", sigma);
        }

        // coefficients based on rotational speed
        fprintf(pFile, "%12s    Ct:%11.5F     Cp:%11.5F    J:%11.5F", "", ct, cp, context.adv * common::pi);
        // coefficients based on forward speed
        fprintf(pFile, "%12s    Tc:%11.5F     Pc:%11.5F  adv:%11.5F", "", tc, pc, context.adv);

        if (context.terse) return;

        // find maximum RE on blade
        double remax = 0.0;
        for (int i = 0; i < context.ii; i++) {
            remax = max(context.re[i], remax);
        }
        double reexp = 1.0;
        if (remax >= 1.0e6) {
            reexp = 6.0;
        } else if (remax >= 1.0e3) {
            reexp = 3.0;
        }

        if (reexp == 1.0) {
            fprintf(pFile, "\n  i  r/R    c/R  beta(deg)   CL      Cd    RE    Mach   effi  effp  na.u/U");
        } else {
            fprintf(pFile, "\n  i  r/R   c/R  beta(deg)  CL     Cd    REx10^%1i Mach   effi  effp  na.u/U", (int)reexp);
        }

        double wa, wt;
        double vw, vaw, utotw, cw, sw, effi;
        double utot, vt, vt_adw, va, va_adw, vd, vd_adw, ci, ci_adv, ci_vt,
               si, si_va, w, w_adv, w_vt, w_va, phi, p_adv, p_vt, p_va;
        double mach, bdeg, xre;
        char schar[1] = {' '};
        for (int i = 0; i < context.ii; i += iadd) {
            // use equivalent prop to define local efficiency
            uvadd(context, context.xi[i], wa, wt);
            vw = context.vwak[i];
            vaw = vw * context.xw[i] / context.adw;
            // Freestream velocity component on equiv prop
            utotw = context.urduct;
            cw = context.xi[i] / context.adv - wt - vw;
            sw = utotw                       + wa - vaw;
            effi = (cw / sw)                 * context.adv / context.xw[i];

            // use real prop to define Mach number
            xoper::cscalc(i, utot, wa, wt,
                          vt, vt_adw,
                          va, va_adw,
                          vd, vd_adw,
                          ci, ci_adv, ci_vt,
                          si, si_va,
                          w, w_adv, w_vt, w_va,
                          phi, p_adv, p_vt, p_va);

            mach = w * context.vel / context.vso;

            bdeg = context.beta[i] * 180. / common::pi;
            xre = context.re[i] / (pow(10., reexp));

            schar[0] = ' ';
            if (context.stall[i]) schar[0] = 'S';

            fprintf(pFile, " %2i%6.3F%7.4F%7.2F%7.3F %1s%7.4F %6.2F %6.3F %6.3F%6.3F%8.3F\n",
                i, context.xi[i], context.ch[i], bdeg, context.cl[i], schar, context.cd[i], xre, mach,
                effi, context.effp[i], context.ubody[i]);
        }
    }

    /**
     * @param context
     * @param xiw
     * @param wa
     * @param wt
     */
    void uvadd(common::context& context, double xiw, double& wa, double& wt) {
        wa = 0;
        wt = 0;

        if (context.nadd <= 1) return;

        double rdim = xiw * context.rad;
        if (rdim >= context.radd[0] and rdim <= context.radd[context.nadd-1]) {
            wa = spline::seval(rdim, context.uadd, context.uaddr, context.radd, context.nadd) / context.vel;
            wt = spline::seval(rdim, context.vadd, context.vaddr, context.radd, context.nadd) / context.vel;
        }
    }
}