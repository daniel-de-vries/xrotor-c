//
// Created by daniel.devries on 11/29/2018.
//

#include <varargs.h>
#include <iostream>
#include <cstdio>
#include <cstdarg>

#include "spline.h"
#include "xaero.h"
#include "userio.h"
#include "xio.h"

namespace xio {

    void load(common::context &context, string fname) {
        context.greek = false;

        if (fname.empty()) {
            userio::asks("Enter filename^", context.fname);
        } else {
            context.fname = fname;
        }
        ifstream ifs(fname);

        string line;

        // File version
        rdline(ifs, line);
        cout << line.substr(16, 21) << endl;

        const char* d4 = "%lf %lf %lf %lf";
        const char* d3 = "%lf %lf %lf";
        const char* d2 = "%lf %lf";
        const char* d1 = "%lf";

        // Case title
        rdline(ifs, context.name);

        if (!read(ifs, line, d4, 4, &context.rho, &context.vso, &context.rmu, &context.alt) or
            !read(ifs, line, d4, 4, &context.rad, &context.vel, &context.adv, &context.rake) or
            !read(ifs, line, d2, 2, &context.xi0, &context.xw0) or
            !read(ifs, line, "%u", 1, &context.naero)) {
            // goto 210
        }

        double xisect,
               a0deg, dclda, clmax, clmin,
               dclda_stall, dcl_stall, cmcon, mcrit,
               cdmin, cldmin, dcddcl2,
               reref, rexp,
               a0;
        for (unsigned n = 0; n < context.naero; n++) {
            if (!read(ifs, line, d1, 1, &xisect) or
                !read(ifs, line, d4, 4, &a0deg, &dclda, &clmax, &clmin) or
                !read(ifs, line, d4, 3, &dclda_stall, &dcl_stall, &cmcon, &mcrit) or
                !read(ifs, line, d3, 3, &cdmin, &cldmin, &dcddcl2) or
                !read(ifs, line, d2, 2, &reref, &rexp)) {
                // goto 210
            } else {
                a0 = a0deg * common::pi / 180.;
                xaero::putaero(context, n, xisect, a0, clmax, clmin,
                               dclda, dclda_stall, dcl_stall,
                               cdmin, cldmin, dcddcl2, cmcon, mcrit, reref, rexp);
            }
        }

        char free, duct;
        if (!read(ifs, line, "%c %c %*c", 3, &free, &duct)) {
            // goto 210
        } else {
            context.free = (free == 'T' or free == 't');
            context.duct = (duct == 'T' or duct == 't');
        }

        cout << endl;
        if (context.free) cout << "Self-deforming wake option set" << endl;
        else              cout << "Rigid wake option set" << endl;
        if (context.duct) cout << "Duct option set" << endl;
        else              cout << "Free-tip option set" << endl;
        cout << endl;

        unsigned iix;
        if (!read(ifs, line, "%u %u", 2, &iix, context.nblds)) {
            // goto 210
        } else {
            double betadeg;
            for (unsigned i = 0; i < iix; i++) {
                if (!read(ifs, line, d4, 4, &context.xi[i], &context.ch[i], &betadeg, &context.ubody[i])) {
                    // goto 210
                } else {
                    context.beta[i] = betadeg * common::pi / 180.;
                    context.beta0[i] = context.beta[i];
                }
            }
        }

        // TODO: Optional duct velocity
        // TODO: Optional slipstream velocities

        // Check for number of analysis sections to use
        if (iix != context.ii) {
            cout << "Read  # input stations = " << iix << endl;
            cout << "Using # blade stations = " << context.ii << endl;
            cout << "Enter # stations or <cr> for " << context.ii << " ";
            cin >> context.ii;
        }

        // spline blade geometry to "old" radial locations
        vec w1(context.xi), w2(context.ch), w3(iix), w4(context.beta), w5(iix), w6(context.ubody), w7(iix);
        spline::spline(w2, w3, w1);
        spline::spline(w4, w5, w1);
        spline::spline(w6, w7, w1);

        // TODO: call xrotor::setx()
        // TODO: call xoper::xwinit()

        // interpolate read-in geometry to generated radial stations
        for (unsigned i = 0; i < context.ii; i++) {
            context.ch[i]       = spline::seval(context.xi[i], w2, w3, w1);
            context.beta[i]     = spline::seval(context.xi[i], w4, w5, w1);
            context.ubody[i]    = spline::seval(context.xi[i], w6, w7, w1);
            context.beta0[i]    = context.beta[i];
        }
        context.iinf = context.ii + context.ii / 2;

        xaero::setiaero(context);
    }

    void rdline(ifstream &ifs, string &line) {
        while (line.find("!#") == string::npos || line.find_first_not_of(" \n") == string::npos) {
            ifs >> line;
            if (ifs.eof()) {
                line = "END";
                return;
            }
        }
    }


    bool read(ifstream &ifs, string &line, const char *fmt, int n, ...) {
        rdline(ifs, line);
        va_list args;
        va_start(args, n);
        bool error = (sscanf(line.c_str(), fmt, args) == n);
        va_end(args);
        return error;
    }
}