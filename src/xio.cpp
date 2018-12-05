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

    void LOAD(common::context &ctxt, string fname) {
        ctxt.GREEK = false;

        if (fname.empty()) {
            userio::ASKS("Enter filename^", ctxt.FNAME);
        } else {
            ctxt.FNAME = fname;
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
        RDLINE(ifs, ctxt.NAME);

        if (!READ(ifs, line, d4, 4, &ctxt.RHO, &ctxt.VSO, &ctxt.RMU, &ctxt.ALT) or
            !READ(ifs, line, d4, 4, &ctxt.RAD, &ctxt.VEL, &ctxt.ADV, &ctxt.RAKE) or
            !READ(ifs, line, d2, 2, &ctxt.XI0, &ctxt.XW0) or
            !READ(ifs, line, "%u", 1, &ctxt.NAERO)) {
            // goto 210
        }

        double xisect,
               a0deg, dclda, clmax, clmin,
               dclda_stall, dcl_stall, cmcon, mcrit,
               cdmin, cldmin, dcddcl2,
               reref, rexp,
               a0;
        for (unsigned n = 0; n < ctxt.NAERO; n++) {
            if (!READ(ifs, line, d1, 1, &xisect) or
                !READ(ifs, line, d4, 4, &a0deg, &dclda, &clmax, &clmin) or
                !READ(ifs, line, d4, 3, &dclda_stall, &dcl_stall, &cmcon, &mcrit) or
                !READ(ifs, line, d3, 3, &cdmin, &cldmin, &dcddcl2) or
                !READ(ifs, line, d2, 2, &reref, &rexp)) {
                // goto 210
            } else {
                a0 = a0deg * common::PI / 180.;
                xaero::PUTAERO(ctxt, n, xisect, a0, clmax, clmin,
                               dclda, dclda_stall, dcl_stall,
                               cdmin, cldmin, dcddcl2, cmcon, mcrit, reref, rexp);
            }
        }

        char free, duct;
        if (!READ(ifs, line, "%c %c %*c", 3, &free, &duct)) {
            // goto 210
        } else {
            ctxt.FREE = (free == 'T' or free == 't');
            ctxt.DUCT = (duct == 'T' or duct == 't');
        }

        cout << endl;
        if (ctxt.FREE) cout << "Self-deforming wake option set" << endl;
        else              cout << "Rigid wake option set" << endl;
        if (ctxt.DUCT) cout << "Duct option set" << endl;
        else              cout << "Free-tip option set" << endl;
        cout << endl;

        unsigned iix;
        if (!READ(ifs, line, "%u %u", 2, &iix, ctxt.NBLDS)) {
            // goto 210
        } else {
            double betadeg;
            for (unsigned i = 0; i < iix; i++) {
                if (!READ(ifs, line, d4, 4, &ctxt.XI[i], &ctxt.CH[i], &betadeg, &ctxt.UBODY[i])) {
                    // goto 210
                } else {
                    ctxt.BETA[i] = betadeg * common::PI / 180.;
                    ctxt.BETA0[i] = ctxt.BETA[i];
                }
            }
        }

        // TODO: Optional DUCT velocity
        // TODO: Optional slipstream velocities

        // Check for number of analysis sections to use
        if (iix != ctxt.II) {
            cout << "Read  # input stations = " << iix << endl;
            cout << "Using # blade stations = " << ctxt.II << endl;
            cout << "Enter # stations or <cr> for " << ctxt.II << " ";
            cin >> ctxt.II;
        }

        // spline blade geometry to "old" radial locations
        vec w1(ctxt.XI), w2(ctxt.CH), w3(iix), w4(ctxt.BETA), w5(iix), w6(ctxt.UBODY), w7(iix);
        spline::SPLINE(w2, w3, w1);
        spline::SPLINE(w4, w5, w1);
        spline::SPLINE(w6, w7, w1);

        // set radial stations for built-in distribution scheme
        xrotor::SETX(ctxt);
        xoper::XWINIT();

        // interpolate read-in geometry to generated radial stations
        for (unsigned i = 0; i < ctxt.II; i++) {
            ctxt.CH[i]       = spline::SEVAL(ctxt.XI[i], w2, w3, w1);
            ctxt.BETA[i]     = spline::SEVAL(ctxt.XI[i], w4, w5, w1);
            ctxt.UBODY[i]    = spline::SEVAL(ctxt.XI[i], w6, w7, w1);
            ctxt.BETA0[i]    = ctxt.BETA[i];
        }
        ctxt.IINF = ctxt.II + ctxt.II / 2;

        xaero::SETIAERO(ctxt);

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