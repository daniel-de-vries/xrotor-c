//
// Created by daniel.devries on 11/29/2018.
//

#include <iostream>
#include <cstdio>
#include <cstdarg>

#include <spline.h>
#include <userio.h>
#include <xaero.h>
#include <xio.h>
#include <xrotor.h>
#include <xrotor.h>

namespace xio {

    void LOAD(common::context &context, string fname) {
        context.GREEK = false;

        if (fname.empty()) {
            userio::ASKS("Enter filename^", context.FNAME);
        } else {
            context.FNAME = fname;
        }
        ifstream ifs(fname);

        string line;

        // File version
        RDLINE(ifs, line);
        cout << line.substr(16, 21) << endl;

        const char* d4 = "%lf %lf %lf %lf";
        const char* d3 = "%lf %lf %lf";
        const char* d2 = "%lf %lf";
        const char* d1 = "%lf";

        // Case title
        RDLINE(ifs, context.NAME);

        if (!READ(ifs, line, d4, 4, &context.RHO, &context.VSO, &context.RMU, &context.ALT) or
            !READ(ifs, line, d4, 4, &context.RAD, &context.VEL, &context.ADV, &context.RAKE) or
            !READ(ifs, line, d2, 2, &context.XI0, &context.XW0) or
            !READ(ifs, line, "%u", 1, &context.NAERO)) {
            // goto 210
        }

        double xisect,
               a0deg, dclda, clmax, clmin,
               dclda_stall, dcl_stall, cmcon, mcrit,
               cdmin, cldmin, dcddcl2,
               reref, rexp,
               a0;
        for (unsigned n = 0; n < context.NAERO; n++) {
            if (!READ(ifs, line, d1, 1, &xisect) or
                !READ(ifs, line, d4, 4, &a0deg, &dclda, &clmax, &clmin) or
                !READ(ifs, line, d4, 3, &dclda_stall, &dcl_stall, &cmcon, &mcrit) or
                !READ(ifs, line, d3, 3, &cdmin, &cldmin, &dcddcl2) or
                !READ(ifs, line, d2, 2, &reref, &rexp)) {
                // goto 210
            } else {
                a0 = a0deg * common::PI / 180.;
                xaero::PUTAERO(context, n, xisect, a0, clmax, clmin,
                               dclda, dclda_stall, dcl_stall,
                               cdmin, cldmin, dcddcl2, cmcon, mcrit, reref, rexp);
            }
        }

        char free, duct;
        if (!READ(ifs, line, "%c %c %*c", 3, &free, &duct)) {
            // goto 210
        } else {
            context.FREE = (free == 'T' or free == 't');
            context.DUCT = (duct == 'T' or duct == 't');
        }

        cout << endl;
        if (context.FREE) cout << "Self-deforming wake option set" << endl;
        else              cout << "Rigid wake option set" << endl;
        if (context.DUCT) cout << "Duct option set" << endl;
        else              cout << "Free-tip option set" << endl;
        cout << endl;

        unsigned iix;
        if (!READ(ifs, line, "%u %u", 2, &iix, context.NBLDS)) {
            // goto 210
        } else {
            double betadeg;
            for (unsigned i = 0; i < iix; i++) {
                if (!READ(ifs, line, d4, 4, &context.XI[i], &context.CH[i], &betadeg, &context.UBODY[i])) {
                    // goto 210
                } else {
                    context.BETA[i] = betadeg * common::PI / 180.;
                    context.BETA0[i] = context.BETA[i];
                }
            }
        }

        // TODO: Optional DUCT velocity
        // TODO: Optional slipstream velocities

        // Check for number of analysis sections to use
        if (iix != context.II) {
            cout << "Read  # input stations = " << iix << endl;
            cout << "Using # blade stations = " << context.II << endl;
            cout << "Enter # stations or <cr> for " << context.II << " ";
            cin >> context.II;
        }

        // spline blade geometry to "old" radial locations
        vec w1(context.XI), w2(context.CH), w3(iix), w4(context.BETA), w5(iix), w6(context.UBODY), w7(iix);
        spline::SPLINE(w2, w3, w1);
        spline::SPLINE(w4, w5, w1);
        spline::SPLINE(w6, w7, w1);

        // set radial stations for built-in distribution scheme
        xrotor::SETX(context);
        xoper::XWINIT();

        // interpolate read-in geometry to generated radial stations
        for (unsigned i = 0; i < context.II; i++) {
            context.CH[i]       = spline::SEVAL(context.XI[i], w2, w3, w1);
            context.BETA[i]     = spline::SEVAL(context.XI[i], w4, w5, w1);
            context.UBODY[i]    = spline::SEVAL(context.XI[i], w6, w7, w1);
            context.BETA0[i]    = context.BETA[i];
        }
        context.IINF = context.II + context.II / 2;

        xaero::SETIAERO(context);

        // XROTOR continues to calculate the operating point here
        // We skip this, since we want to use OPER to calculate an off-design operating point
    }

    void RDLINE(ifstream &ifs, string &line) {
        while (line.find("!#") == string::npos || line.find_first_not_of(" \n") == string::npos) {
            ifs >> line;
            if (ifs.eof()) {
                line = "END";
                return;
            }
        }
    }


    bool READ(ifstream &ifs, string &line, const char *fmt, int n, ...) {
        RDLINE(ifs, line);
        va_list args;
        va_start(args, n);
        bool error = (sscanf(line.c_str(), fmt, args) == n);
        va_end(args);
        return error;
    }
}