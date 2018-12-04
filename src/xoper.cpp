//
// Created by daniel.devries on 12/3/2018.
//

#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

#include <vortex.h>
#include <userio.h>
#include <xaero.h>
#include <xoper.h>
#include <xrotor.h>
#include <xutils.h>

using common::fprintf;

namespace xoper {

    /**
     * Run rotor at arbitrary operating points.
     *
     * @param context
     */
    void OPER(common::context &context) {
        double plfac1 = 0.7;
        double plfac2 = 0.8;
        double plfacd = 0.6;
        double xorg = 0.15;
        double yorg = 0.10;

        context.GREEK = false;

        bool error;
        int ninput;
        vector<int> iinput(20);
        vector<double> rinput(20);

        bool showHelp = false;
        while (true) {
            if (showHelp) {
                printf(
//                 "\n   ADVA r   Prescribe advance ratio"
                 "\n   RPM  r   Prescribe rpm"
                 "\n   THRU r   Prescribe thrust"
//                 "\n   TORQ r   Prescribe torque"
//                 "\n   POWE r   Prescribe power"
//                 "\n"
//                 "\n   ASEQ rrr Calculate case sequence of advance ratios"
//                 "\n   RSEQ rrr Calculate case sequence of rpms"
//                 "\n   BSEQ rrr Calculate case sequence of blade angles"
//                 "\n   VSEQ rrr Calculate case sequence of speeds at fixed pitch"
//                 "\n   CLRC     Clear case accumulator"
//                 "\n   ADDC     Add current point point to case accumulator"
//                 "\n   CPUT f   Write current case accumulator to file"
//                 "\n   CGET f   Read cases from file"
//                 "\n   CASE i   Select case"
//                 "\n"
//                 "\n   ATMO r   Set fluid properties from standard atmosphere"
//                 "\n   VELO r   Set or change flight speed"
//                 "\n   ANGL r   Change blade pitch angle"
//                 "\n   PVAR f   Enter and use engine rpm/power line"
                 "\n"
                 "\n   FORM     Select slipstream and velocity formulation"
                 "\n"
//                 "\n   NAME s   Set or change case name"
                 "\n   WRIT f   Write current operating point to disk file"
                 "\n   DISP     Display current operating state"
//                 "\n   INIT     Initialize next analysis case"
//                 "\n   REIN     Re-initialize prop to known operating state"
//                 "\n   TERS     Toggle between TERSE and verbose OUTPUT"
                 "\n   ITER i   Change max number of Newton iterations"
                 "\n   N    i   Change number of radial points");
                showHelp = false;
            }

            string command, comarg;
            userio::ASKC(".OPER", command, comarg);

            for (int i = 0; i < 20; i++) {
                iinput[i] = 0;
                rinput[i] = 0.0;
            }
            ninput = 0;
            userio::GETINT(comarg, iinput, ninput, error);
            ninput = 0;
            userio::GETFLT(comarg, rinput, ninput, error);

                 if (command == "    ") return;
            else if (command == "?   ") {showHelp = true; continue;}
            else if (command == "FORM") commands::FORM(context);    // goto 2
//            else if (command == "TERS") commands::TERS(context);    // goto 4
            else if (command == "DISP") commands::DISP(context);    // goto 10
//            else if (command == "NAME") commands::NAME(context);    // goto 15
            else if (command == "WRIT") commands::WRIT(context, comarg);    // goto 20
//            else if (command == "DUCT") commands::DUCT(context);    // goto 22
//            else if (command == "VRAT") commands::VRAT(context);    // goto 24
//            else if (command == "ATMO") commands::ATMO(context);    // goto 35
//            else if (command == "VELO") commands::VELO(context);    // goto 38
//            else if (command == "ANGL") commands::ANGL(context);    // goto 40
//            else if (command == "ADVA") commands::ADVA(context);    // goto 42
            else if (command == "RPM ") commands::RPM (context, ninput, rinput);    // goto 45
            else if (command == "THRU") commands::THRU(context, ninput, rinput);    // goto 50
//            else if (command == "TORQ") commands::TORQ(context);    // goto 60
//            else if (command == "POWE") commands::POWE(context);    // goto 70
//            else if (command == "ASEQ") commands::ASEQ(context);    // goto 81
//            else if (command == "RSEQ") commands::RSEQ(context);    // goto 82
//            else if (command == "BSEQ") commands::BSEQ(context);    // goto 83
//            else if (command == "VSEQ") commands::VSEQ(context);    // goto 84
//            else if (command == "CLRC") commands::CLRC(context);    // goto 90
//            else if (command == "ADDC") commands::ADDC(context);    // goto 92
//            else if (command == "CPUT") commands::CPUT(context);    // goto 94
//            else if (command == "CGET") commands::CGET(context);    // goto 96
//            else if (command == "CASE") commands::CASE(context);    // goto 97
//            else if (command == "LIST") commands::LIST(context);    // goto 98
//            else if (command == "N   ") commands::N   (context);    // goto 72
            else if (command == "ITER") commands::ITER(context, ninput, iinput);    // goto 75
//            else if (command == "INIT") commands::INIT(context);    // goto 76
//            else if (command == "REIN") commands::REIN(context);    // goto 78
            else {
                printf(" %4s command not recognized.\n  Type \"?\" for list, <Return> to exit menu.", command.c_str());
            }
        }
    }

    namespace commands {
        /**
         * Select options for slipstream and velocity calculation.
         * @param context
         */
        void FROM(common::context& context) {
            string command, comarg;
            while (true) {
                printf(
                    "\n Select options for calculation of slipstream velocities"
                    "\n   GRAD     use Graded Momentum       Formulation "
                    "\n   POT      use Potential (Goldstein) Formulation "
                    "\n   VRTX     use discrete Vortex Wake  Formulation "
                    "\n   WAKE     Toggle between rigid and self-deforming wake");

                userio::ASKC(".FORM", command, comarg);
                     if (command == "GRAD") {context.VRTX = false; context.FAST = true;}
                else if (command == "POT ") {context.VRTX = false; context.FAST = false;}
                else if (command == "VRTX") {context.VRTX = true;}
                else if (command == "WAKE") {context.FREE = !context.FREE;}
                else if (command == "    ") {return;}

                if (context.VRTX) printf("Discrete Vortex Formulation selected\n");
                else if (context.FAST) printf("Graded Momentum Formulation selected\n");
                else printf("Potential Formulation selected\n");

                if (context.FREE) printf("Self-deforming wake selected");
                else printf("Rigid wake selected");
            }
        }

        /**
         * Display current prop operating point data.
         * @param context
         */
        void DISP(common::context& context) {
            xrotor::OUTPUT(context, cout);
        }

        /**
         * Write current prop operating point data to file
         * @param context
         * @param fName         filename
         */
        void WRIT(common::context& context, string fName) {
            if (fName[0] != ' ') context.SAVIL = fName;
            ofstream ofs;
            xrotor::OPFILE(ofs, context.SAVIL);
            xrotor::OUTPUT(context, ofs);
            ofs.close();
        }

        /**
         * Specify RPM and solve.
         * @param context
         * @param ninput        number of input arguments passed to RPM
         * @param rinput        list of input arguments passed to RPM
         */
        void RPM(common::context& context, int ninput, vector<double> rinput) {
            double rpm;
            if (ninput >= 1) {
                rpm = rinput[0];
            } else {
                rpm = context.VEL / (context.RAD * context.ADV * common::PI / 30.);
                userio::ASKR("rpm               ", rpm);
            }
            context.ADV = context.VEL / (context.RAD * rpm * common::PI / 30.);
            context.CONV = false;
            APER(context, 4, 2, context.LOPRINI);
            if (context.CONV) xrotor::OUTPUT(context, cout);
        }

        /**
         * Specify thrust and solve.
         * @param context
         * @param ninput        number of input arguments passed to THRU
         * @param rinput        list of input arguments passed to THRU
         */
        void THRU(common::context& context, int ninput, vector<double> rinput) {
            if (ninput >= 1) {
                context.TSPEC = rinput[0];
            } else {
                context.TSPEC = context.TTOT * context.RHO * pow(context.VEL, 2) * pow(context.RAD, 2);
                userio::ASKR("thrust (N)        ", context.TSPEC);
            }
            double rpm = context.VEL / (context.RAD * context.ADV * common::PI / 30.);
            printf("\n Current rpm:%9.2f", rpm);
            string ans, ansarg;
            while (ans != "R   " and ans != "P   " ) {
                userio::ASKC("fix Pitch / fix Rpm ( P/R )?", ans, ansarg);
            }

            context.CONV = false;
            double bsav = context.BETA[context.II-1];
            if (ans == "P   ") APER(context, 1, 2, context.LOPRINI);
            else {
                userio::ASKR("rpm:", rpm);
                context.ADV = context.VEL / (context.RAD * rpm * common::PI / 30.);
                APER(context, 1, 1, context.LOPRINI);
            }

            if (context.CONV) {
                xrotor::OUTPUT(context, cout);
            }
            // Check for valid blade angle change
            if (ans != "P   ") {
                if (context.CONV) {
                    // convergence was achieved: show blade angle change incurred
                    printf(" Blade angle changed %7.3f degrees", context.DBETA * 180. / common::PI);
                } else {
                    // convergence failed: restore clobbered blade angles
                    for (int i = 0; i < context.II; i++) {
                        context.BETA[i] -= context.DBETA;
                        context.BETA0[i] -= context.DBETA;
                    }
                }
            }
        }

        /**
         * Set max number or iterations for nonlinear solution
         * @param context
         * @param ninput        number of input arguments passed to ITER
         * @param iinput        list of input arguments passed to ITER
         */
        void ITER(common::context& context, int ninput, vector<int> iinput) {
            if (ninput >= 1) {
                context.NITERA = iinput[0];
            } else {
                userio::ASKI("Max number of iterations", context.NITERA);
            }
        }
    }

    /**
     * Set reasonable initial circulation and converge arbitrary operating point.
     *
     * @param context
     * @param ispec, icon   @see[xoper::apiter()]
     * @param linit         flag for initialization of rotor conditions
     */
    void APER(common::context &context, unsigned short ispec, unsigned short icon, bool linit) {
        // Initialize circulations if requested
        if (linit) APINIT(context);

        APITER(context, ispec, icon);

        if (!context.CONV) {
            cout << endl;
            cout << "Iteration limit exceeded" << endl;
            cout << "Gres Fres Ares = " << context.GRESMX << " " << context.FRESMX << " " << context.ARESMX << endl;
        }
    }

    /**
     * Set reasonable initial circulation.
     *
     * Initial circulations are set w/o induced effects.
     * An iteration is done using the self-induced velocity
     * from graded momentum theory to converge an approximate
     * wake advance ratio.
     *
     * @param context
     */
    void APINIT(common::context &context) {
        const unsigned niterg = 10;

        double blds = (float) context.NBLDS;
        context.DBETA = 0.0;

        double uduct     = 0.0;
        double vaduct_va;
        if (context.DUCT) {
            uduct     = context.URDUCT - 1;
        }
        context.ADW = context.ADV * (1.0 + uduct);

        // ==========================================================
        // Initialize section circulation neglecting induced velocity
        double tsum = 0.;
        double utot, wa, wt, si, ci, wsq, w, phi, alfa, rey;
        double cl_al, cl_w, clmax, clmin, dclstall, cd_alf, cd_w, cd_rey, cm_al, cm_w;
        for (int i = 0; i < context.II; i++) {
            utot = context.URDUCT + context.UBODY[i];
            xrotor::UVADD(context, context.XI[i], wa, wt);

            si = utot                        + wa;
            ci = context.XI[i] / context.ADV - wt;

            wsq = ci*ci + si*si;
            w = sqrt(wsq);
            phi = atan2(si, ci);

            alfa = context.BETA[i] - phi;
            rey = context.CH[i] * abs(w) * context.RHO * context.VEL * context.RAD / context.RMU;
            xaero::GETCLCDCM(context, i, alfa, w, rey,
                             context.CL[i], cl_al, cl_w,
                             clmax, clmin, dclstall, context.STALL[i],
                             context.CD[i], cd_alf, cd_w, cd_rey,
                             context.CM[i], cm_al, cm_w);

            context.GAM[i] = 0.5 * context.CL[i] * w * context.CH[i];
            tsum += blds * context.GAM[i] * ci * context.DXI[i];
        }

        // use momentum theory estimate of axial velocity to set wake ADV. ratio
        double vhsq = 0.5 * tsum / common::PI;
        vhsq = max(vhsq, -0.25);
        context.ADW = context.ADV * 0.5 * (1.0 + sqrt(1.0 + 4.0 * vhsq));

        // recalculate Vtan using new GAM values
        VCALC(context);
        // ==========================================================

        // ===============================================================
        // Refine the initial guess with a graded-momentum theory estimate
        // Use momentum theory to estimate axial induced velocity to drive
        // equation for wake advance ratio
        double t_adw, dclmax, rlxmin, vt, vt_gam, vt_adw, va, va_gam, va_adw;
        double ci_adv, ci_vt, si_va, w_adv, w_vt, w_va, p_adv, p_vt, p_va, al_vt, al_va;
        double rez, z_cl, z_w, z_g, z_al, z_vt, z_va, z_adw, delg, dcl;
        double rlx, g_adw, cosr, t_g, t_vt, vhsq_t;
        for (int iterg = 0; iterg < niterg; iterg++) {
            GRADMO(context.II, context.NBLDS, context.DUCT, context.RAKE,
                   context.XI, context.XV, context.GAM, context.ADW, context.VIND_GAM, context.VIND_ADW);

            tsum = 0;
            t_adw = 0;

            dclmax = 0;
            rlxmin = 0;

            for (int i = 0; i < context.II; i++) {
                // Redefine VT and VA to diagonal self-influences
                vt     = context.VIND_GAM[2][i][i] * context.GAM[i];
                vt_gam = context.VIND_GAM[2][i][i];
                vt_adw = context.VIND_ADW[2][i];

                va     = context.VIND_GAM[0][i][i] * context.GAM[i];
                va_gam = context.VIND_GAM[0][i][i];
                va_adw = context.VIND_ADW[0][i];

                // include duct effect on freestream and induced axial velocity
                uduct = 0.0;
                vaduct_va      = 1.0;
                if (context.DUCT) {
                    uduct = context.URDUCT - 1.0;
                    vaduct_va = 2.0 * context.URDUCT;
                }

                utot = 1.0 + uduct + context.UBODY[i];
                xrotor::UVADD(context, context.XI[i], wa, wt);

                ci     =  context.XI[i] / context.ADV - wt - vt;
                ci_adv = -context.XI[i] / pow(context.ADV, 2);
                ci_vt  = - 1.0;

                si     = utot + wa + va * vaduct_va;
                si_va  = vaduct_va;

                wsq = ci*ci + si*si;
                w = sqrt(wsq);
                w_adv = (ci * ci_adv              ) / w;
                w_vt  = (ci * ci_vt               ) / w;
                w_va  = (              si * si_va ) / w;

                phi = atan2(si, ci);
                p_adv = (            - si * ci_adv) / wsq;
                p_vt  = (            - si * ci_vt ) / wsq;
                p_va  = (ci * si_va               ) / wsq;

                alfa  = context.BETA[i] - phi;
                al_vt =                 - p_vt;
                al_va =                 - p_va;

                rey = context.CH[i] * abs(w) * context.RHO * context.VEL * context.RAD / context.RMU;
                xaero::GETCLCDCM(context, i, alfa, w, rey,
                                 context.CL[i], cl_al, cl_w,
                                 clmax, clmin, dclstall, context.STALL[i],
                                 context.CD[i], cd_alf, cd_w, cd_rey,
                                 context.CM[i], cm_al, cm_w);

                // Res( CL( AL W ) , W , GAM )
                rez   = context.CH[i] * context.CL[i] * w - 2.0 * context.GAM[i];
                z_cl  = context.CH[i]                 * w;
                z_w   = context.CH[i] * context.CL[i];
                z_g   =                                   - 2.0;

                // Res( AL( VT ADW ) , W( VT ADW ) , GAM )
                z_al  = z_cl * cl_al;
                z_w   = z_cl * cl_w   + z_w;

                // Res( VT(GAM ADW) , ADW , GAM )
                z_vt  = z_w * w_vt    + z_al * al_vt;
                z_va  = z_w * w_va    + z_al * al_va;

                // Res( ADW , GAM )
                z_adw = z_vt * vt_adw + z_va * va_adw;
                z_g   = z_vt * vt_gam + z_va * va_gam + z_g;

                delg = -rez / z_g;
                dcl = 2.0 * delg / (context.CH[i] * w);

                // Apply limiter to GAM update based on CL change
                rlx = 1.0;
                if (rlx * abs(dcl) > 0.2) rlx = min(rlx, 0.2 / abs(dcl));

                if (abs(dcl) > abs(dclmax)) dclmax = dcl;
                if (abs(rlx) < rlxmin)      rlxmin = rlx;

                context.GAM[i] += rlx * delg;
                // dREZ = Z_G*dG + Z_ADW*dADW = 0
                g_adw = -z_adw / z_g;

                // Forces for raked blade corrected for COS of rake angle
                cosr = 1.0;

                tsum  += blds * context.GAM[i] * ci      * context.DXI[i] * cosr;
                t_g    = blds                  * ci      * context.DXI[i] * cosr;
                t_vt   = blds * context.GAM[i] * ci * vt * context.DXI[i] * cosr;
                t_adw += (t_g + t_vt * vt_gam) * g_adw
                        +       t_vt * vt_adw;
            }

            // Momentum theory estimate of induced axial velocity
            vhsq   = 0.5 * tsum / common::PI;
            vhsq   = max(vhsq, -0.2499);
            vhsq_t = 0.5        / common::PI;

            rez    = context.ADW - context.ADV  * 0.5 * (1.0            + sqrt(1.0 + 4.0 * vhsq));
            z_adw  = 1.0         - context.ADV / sqrt(1.0 + 4.0 * vhsq) * vhsq_t * t_adw;
            if (z_adw == 0) cout << "APRINT Z_ADW " << z_adw << endl;

            context.DADW = -rez / z_adw;
            context.DADW = min(context.DADW, 10.0 * context.ADW);
            context.DADW = max(context.DADW, -0.9 * context.ADW);
            context.ADW += context.DADW;

            if (rlxmin < 0.2) FILTER(context.GAM, 0.2*context.II, context.II);

            if (abs(dclmax) < 0.001) return;
        }
    }

    /**
     * Converge arbitrary performance operating point.
     *
     * @param context
     * @param ispec         controls the quantity used as target quantity:
     *                          - ispec = 1 : drive the thrust to TSPEC
     *                          - ispec = 2 : drive the torque to QSPEC
     *                          - ispec = 3 : drive the power  to PSPEC
     *                          - ispec = 4 : fix advance ratio to current value
     *                          - ispec = 5 : drive to power specified by RPM (engine power-RPM line)
     * @param icon          controls the constrained quantity:
     *                          - icon = true  : advance ratio(rpm) fixed
     *                          - icon = false : blade pitch fixed
     */
    void APITER(common::context &context, unsigned short ispec, bool icon) {
        double clmax[common::IX], clmin[common::IX], dclstall[common::IX];

        // convergence tolerance
        const double eps = 1.0e-7;

        int k1 = context.II + 1;
        int k2 = context.II + 2;
        int k3 = context.II + 3;
        cout << "2000";

        double utot, wa, wt, vt75, vt_adw, va75, va_adw, vd75, vd_adw, ci75, ci_adv, ci_vt, si75, si_va;
        double w75, w_adv, w_vt, w_va, phi75, p_adv, p_vt, p_va, advfact, z_tw, z_pw;
        double vt, va, vd, ci, si, w, phi;
        double alfa, al_dbe, al_p, rey;
        double cl_al, cl_w, cd_alf, cd_w, cd_rey, cm_al, cm_w;
        double z_cl, z_w, z_gi, z_vt, z_va, z_adv, z_dbe;
        double t_spec, q_spec, p_spec;
        for (int iter = 0; iter < max(context.NITERA, 1); iter++) {
            // if wake advance ratio changed, recalculate Vtan influence coefficients
            if (context.FREE or iter == 0) {
                if (context.FAST) {
                    GRADMO(context.II, context.NBLDS, context.DUCT, context.RAKE,
                           context.XI, context.XV, context.GAM, context.ADW, context.VIND_GAM, context.VIND_ADW);

                    context.IWTYP = 1;
                } else if (not context.VRTX) {
                    HELICO(context.II, context.NBLDS, context.DUCT, context.RAKE,
                           context.XI, context.XV, context.GAM, context.ADW, context.VIND_GAM, context.VIND_ADW);

                    context.IWTYP = 2;
                } else {
                    vortex::VRTXC0(context.II, context.NBLDS, context.DUCT, context.RAKE,
                                   context.XI, context.XV, context.GAM, context.ADW, context.VIND_GAM, context.VIND_ADW);

                    context.IWTYP = 3;
                }
            }

            // recalculate Vtan
            VCALC(context);

            // recalculate wake radius array and Vwak
            SETXW(context);

            // recalculate thrust, power, and sensitivities for current solution
            TPQ(context, 1);

            // initialize max residuals
            context.GRESMX = 0.;
            context.FRESMX = 0.;
            context.ARESMX = 0.;

            for (int j = 0; j < k3; j++) {
                context.Q[k2][j] = 0.;
            }

            // The wake advance ratio equation is only approximate, normally the
            // tangential induced velocity is ignored (inconsistent with a rigid
            // wake with constant wake advance ratio).  This calculates a factor
            // to compensate for the Vt term at one (representative) radial station
            int i;
            for (i = 0; i < context.II; i++) {
                if (context.XI[i] > 0.75) break;
            }
            int i75 = i;
            CSCALC(context,
                    i75, utot, wa, wt,
                    vt75, vt_adw,
                    va75, va_adw,
                    vd75, vd_adw,
                    ci75, ci_adv, ci_vt,
                    si75,               si_va,
                     w75,  w_adv,  w_vt, w_va,
                   phi75,  p_adv,  p_vt, p_va);
            // Factor for OMEG*R-VT correction to wake advance ratio
            // advfact = 1.0 / (1.0 - context.ADV * vt75 / context.XI[i75]);
            // Set to 1.0 for now... HHY
            advfact = 1.0;

            if (context.FREE) {
                // Set up equation to converge wake advance ratio based on
                // average axial velocity consistent with basic momentum theory

                // Use "equivalent" prop thrust and power
                context.DQ[k2] =  context.ADWFCTR * context.ADW * context.TWAK / context.PWAK - context.ADV * advfact;
                z_tw           =  context.ADWFCTR * context.ADW                / context.PWAK;
                z_pw           = -context.ADWFCTR * context.ADW * context.TWAK / pow(context.PWAK, 2);
                for (int j = 0; j < context.II; j++) {
                    context.Q[k2][j] = z_tw * context.TW_GAM[j] + z_pw * context.PW_GAM[j];
                }
                context.Q[k2][k2]    = z_tw * context.TW_ADV + z_pw * context.PW_ADV - advfact;
                context.Q[k2][k2]    = z_tw * context.TW_ADV + z_pw * context.PW_ADV - advfact * context.TWAK / context.PWAK;
                context.ARESMX = max(context.ARESMX, abs(context.DQ[k2] / context.ADV));
            } else {
                // specify zero change of wake advance ratios
                context.DQ[k2] = 0.;
                context.Q[k2][k2] = 1.0;
            }

            // go over stations, enforcing Gamma-CL relation at real prop
            for (i = 0; i < context.II; i++) {
                CSCALC(context,
                        i, utot, wa, wt,
                       vt, vt_adw,
                       va, va_adw,
                       vd, vd_adw,
                       ci, ci_adv, ci_vt,
                       si,                si_va,
                        w,  w_adv,  w_vt,  w_va,
                       phi, p_adv,  p_vt,  p_va);

                alfa   = context.BETA[i] - phi;
                al_dbe =  1.0;
                al_p   = -1.0;

                rey = context.CH[i] * abs(w) * context.RHO * context.VEL * context.RAD / context.RMU;
                xaero::GETCLCDCM(context,
                                 i, alfa, w, rey,
                                 context.CL[i], cl_al, cl_w,
                                 clmax[i], clmin[i], dclstall[i], context.STALL[i],
                                 context.CD[i], cd_alf, cd_w, cd_rey,
                                 context.CM[i], cm_al, cm_w);

                // Enforce local Gamma-CL relation
                context.DQ[i] = context.CH[i] * context.CL[i] * w - 2.0 * context.GAM[i];   // Residual
                z_cl          = context.CH[i]                 * w;
                z_w           = context.CH[i] * context.CL[i];

                z_gi  =                     - 2.0;
                z_vt  = z_cl * (cl_al * al_p * p_vt  + cl_w * w_vt )  + z_w * w_vt;
                z_va  = z_cl * (cl_al * al_p * p_va  + cl_w * w_va )  + z_w * w_va;
                z_adv = z_cl * (cl_al * al_p * p_adv + cl_w * w_adv)  + z_w * w_adv;
                z_dbe = z_cl * (cl_al * al_dbe                     );

                for (int j = 0; j < context.II; j++) {
                    context.Q[i][j] = z_vt * context.VIND_GAM[2][i][j]
                                    + z_va * context.VIND_GAM[0][i][j];                     // dRes/dGamj
                }
                context.Q[i][i]  = context.Q[i][i] + z_gi;                                  // dRes/dGami
                context.Q[i][k1] = z_adv;                                                   // dRes/dAdv
                context.Q[i][k2] =                   z_vt * vt_adw + z_va * va_adw;         // dRes/dAdw
                context.Q[i][k3] = z_dbe;                                                   // dRes/dBeta

                context.GRESMX = max(context.GRESMX, abs(context.DQ[i] / (0.1 * w)));
            }

            // equivalent prop will be used to define inviscid thrust
            switch (ispec) {
                case 1: // drive thrust to specified value
                    t_spec = context.TSPEC / (context.RHO * pow(context.VEL, 2) * pow(context.RAD, 2));
                    context.DQ[k1] = context.TWAK + context.TVIS - t_spec;
                    for (int j = 0; j < context.II; j++) {
                        context.Q[k1][j] = context.TW_GAM[j] + context.TV_GAM[j];
                    }
                    context.Q[k1][k1] = context.TW_ADV + context.TV_ADV;
                    context.Q[k1][k2] = context.TW_ADW + context.TV_ADW;
                    context.Q[k1][k3] =                  context.TV_DBE;

                    context.FRESMX  = max(context.FRESMX, abs(context.DQ[k1]));
                case 2: // drive torque (= PTOT*ADV) to specified value
                    // TODO: implement
                    break;
                case 3: // drive power to specified value
                    // TODO: implement
                    break;
                case 4: // fix advance ratio
                    // TODO: implement
                    break;
                case 5: // drive power to value given by RPM
                    // TODO: implement
                    break;
                default: break;
            }

            // Constraint conditions
            context.DQ[k3] = 0.;
            for (int j = 0; j < k3; j++) {
                context.Q[k3][j] = 0.;
            }
            if (icon) context.Q[k3][k1] = 1.0;      // advance ratio(rpm) fixed
            else      context.Q[k3][k3] = 1.0;      // blade pitch fixed

            // solve linearized Newton system
            xutils::GAUSS(k3, context.Q, &context.DQ, 1);
        }
    }

    /**
     * Calculate cartesian induced velocities.
     * @param context
     */
    void VCALC(common::context &context) {
        double vsum[3] = {0, 0, 0};
        for (int i = 0; i < context.II; i++) {
            for (int j = 0; j < context.II; j++) {
                for (int k = 0; k < 3; k++) {
                    vsum[k] += context.VIND_GAM[k][i][j] * context.GAM[i];
                }
            }
            for (int k = 0; k < 3; k++) {
                context.VIND[k][i] = vsum[k];
                vsum[k] = 0;
            }
        }
    }

    /**
     * Calculate "Graded Momentum" Gamma-swirl influence coefficients.
     *
     * Inputs:
     * @tparam imax         array dimension
     * @param ii            number of radial points on blade (circulation stations)
     * @param nblds         number of blades
     * @param lduct         true for duct outer BC
     * @param rake
     * @param xi            radial coordinate array
     * @param xv
     * @param gam           circulation array
     * @param adw           wake advance ratio  V/wR
     *
     * Outputs:
     * @param vind_gam      sensitivity of velocity at i to circulation at j
     * @param vind_adw      sensitivity of velocity at i to wake advance ratio
     *
     * Where vind_xxx(1,i,j) is the axial component
     *       vind_xxx(3,i,j) is the swirl component
     */
    template <const int imax>
    void GRADMO(const int &ii, const int &nblds, const bool &lduct, const double &rake,
                const double xi[imax], const double xv[imax], const double gam[imax], const double &adw,
                double vind_gam[3][imax][imax], double vind_adw[3][imax]) {
        auto blds = (double) nblds;

        double xi0 = xv[0];
        double xitip = xv[ii+1];

        if (lduct) {
            // Circulation defines mean swirl at blade
            // use simple mean swirl to get swirl at blade
            for (int i = 0; i < ii; i++) {
                for (int j = 0; j < ii; j++) {
                    for (int k = 0; k < 3; k++) {
                        vind_gam[k][i][j] = 0.;
                    }
                }
                vind_gam[2][i][i] =  blds / (4.0 * common::PI * xi[i]);
                vind_adw[2][i]    =  0.0;
                vind_adw[1][i]    =  0.0;
                vind_gam[0][i][i] =  vind_gam[2][i][i] * xi[i] / adw;
                vind_adw[0][i]    = -vind_gam[0][i][i] * gam[i] / adw;
            }
        } else {
            // Circulation defines mean swirl at blade
            // Free-tip treatment incorporates Prandtl's averaging factor F
            double sfac = sqrt(1.0 + 1.0 / pow(adw, 2));
            double sf_adw = 0.5 / sfac * (-2.0 / pow(adw, 3));

            double arg, ek, ek_adw, fk, fk_adw, f, f_adw;

            for (int i = 0; i < ii; i++) {
                for (int k = 0; k < 3; k++) {
                    for (int j = 0; j < ii; j++) {
                        vind_gam[k][i][j] = 0.;
                    }
                    vind_adw[k][i] = 0.;
                }

                arg = min(20.0, 0.5 * blds * (1.0 - xi[i] / xitip) * sfac);
                ek = exp(-arg);
                ek_adw = -ek * 0.5 * blds * (1.0 - xi[i] / xitip) * sf_adw;
                fk = sqrt(1.0 - ek*ek);
                fk_adw = 0.5 / fk * (-2.0 * ek * ek_adw);
                f = atan2(fk, ek) * 2.0 / common::PI;
                f_adw = (ek * fk_adw - fk * ek_adw) / (ek*ek + fk*fk) * 2.0 / common::PI;

                vind_gam[2][i][i] = blds / (4.0 * common::PI * f * xi[i]);
                vind_adw[2][i]    = vind_gam[2][i][i] * gam[i] * (-f_adw / f);
                vind_gam[0][i][i] = vind_gam[2][i][i] * xi[i]  / adw;
                vind_adw[0][i]    = vind_adw[2][i]    * xi[i]  / adw
                                   -vind_gam[0][i][i] * gam[i] / adw;
            }
        }
    }
}