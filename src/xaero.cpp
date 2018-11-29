//
// Created by daniel.devries on 11/29/2018.
//

#include <algorithm>    // std::max, std::min
#include <cmath>        // std::abs, std::pow, std::sqrt, std::exp, std::log
#include <iostream>     // std::cout, std::endl
#include "xaero.h"
using std::max;
using std::min;
using std::cout;
using std::endl;
using std::abs;
using std::pow;
using std::sqrt;
using std::exp;
using std::log;

namespace xaero {

    /**
     * Set up indices referring to aero section for each radial section.
     *
     * @param context
     */
    void setiaero(common::context &context) {
        context.iaero.resize(context.ii);
        for (unsigned long i = 0; i < context.ii; i++) {
            context.iaero[i] = 0;
            for (unsigned n = 0; n < context.naero; n++) {
                if (context.xiaero[n] <= context.xi[i])
                    context.iaero[i] = n;
            }
        }
    }

    /**
     * Put aero data into stored section array at index n.
     *
     * @param context
     * @param n             section index
     * @param xisect        r/R of the section
     * @param a0            angle of zero lift
     * @param clmax         maximum lift coefficient
     * @param clmin         minimum lift coefficient
     * @param dclda         incompressible 2D lift curve slope
     * @param dclda_stall   2D lift curve slope at stall
     * @param dcl_stall     lift coefficient increment, onset to full stall
     * @param cdmin         minimum drag coefficient
     * @param clcdmin       lift at minimum drag coefficient
     * @param dcddcl2       parabolic drag parameter, d(C_d)/d(C_l^2)
     * @param cmcon         incompressible 2D pitching moment
     * @param mcrit         critical Mach number
     * @param reref         reference Reynold's number
     * @param rexp          Reynold's number exponent (Cd ~ Re^rexp)
     */
    void putaero(common::context &context,
                 unsigned n, double xisect, double a0, double clmax, double clmin,
                 double dclda, double dclda_stall, double dcl_stall,
                 double cdmin, double clcdmin, double dcddcl2,
                 double cmcon, double mcrit, double reref, double rexp) {
        if (n > context.aerodata.size())
            context.aerodata.resize(n);
        if (n > context.xiaero.size())
            context.xiaero.resize(n);

        context.aerodata[n][0 ] = a0;
        context.aerodata[n][1 ] = clmax;
        context.aerodata[n][2 ] = clmin;
        context.aerodata[n][3 ] = dclda;
        context.aerodata[n][4 ] = dclda_stall;
        context.aerodata[n][5 ] = dcl_stall;
        context.aerodata[n][6 ] = cdmin;
        context.aerodata[n][7 ] = clcdmin;
        context.aerodata[n][8 ] = dcddcl2;
        context.aerodata[n][9 ] = cmcon;
        context.aerodata[n][10] = reref;
        context.aerodata[n][11] = rexp;
        context.aerodata[n][12] = mcrit;
        context.xiaero[n]       = xisect;
    }

    /**
     * Interpolate C_l(alpha), C_d(alpha), and C_m(alpha) at blade station with index is.
     *
     * @param context
     * @param is                        blade station index
     * @param alf                       angle of attack, alpha
     * @param w                         ?
     * @param rey                       Reynolds number, Re
     * @param clift                     lift coefficient, C_l
     * @param cl_alf                    d(C_l)/d(alpha)
     * @param cl_w                      d(C_l)/d(w)
     * @param clmax, clmin, dcl_stall   @see[xaero::putaero()]
     * @param stallf                    true if stalled
     * @param cdrag                     drag coefficient, C_d
     * @param cd_alf                    d(C_d)/d(alpha)
     * @param cd_w                      d(C_d)/d(w)
     * @param cd_rey                    d(C_d)/d(Re)
     * @param cmom                      moment coefficient, C_m
     * @param cm_al                     d(C_m)/d(alpha)
     * @param cm_w                      d(C_m)/d(w)
     */
    void getclcdcm(common::context &context,
                   unsigned is, double alf, double w, double rey,
                   double &clift, double &cl_alf, double &cl_w,
                   double &clmax, double &clmin, double &dcl_stall, bool &stallf,
                   double &cdrag, double &cd_alf, double &cd_w, double &cd_rey,
                   double &cmom, double &cm_al, double &cm_w) {
        // Check for installed aero data section index
        unsigned n = context.iaero[is];
        if (n < 0 || n >= context.naero) {
            bool error = false;

            if (context.naero > 0) {
                // Find lower index of aero data sections XIAERO(N) bounding XI(IS)
                for (n = 0; n < context.naero; n++) {
                    if (context.xiaero[n] <= context.xi[is]) {
                        context.iaero[is] = n;
                    } else {
                        error = true;
                        break;
                    }
                }

                cout << "Aero section not found for station " << context.xi[is] << endl;
            }

            if (error) {
                n = 0;
                context.iaero[is] = n;
            }
        }

        // Get section aero data from stored section array
        double a0           = context.aerodata[n][0 ];
        clmax               = context.aerodata[n][1 ];
        clmin               = context.aerodata[n][2 ];
        double dclda        = context.aerodata[n][3 ];
        double dclda_stall  = context.aerodata[n][4 ];
        dcl_stall           = context.aerodata[n][5 ];
        double cdmin        = context.aerodata[n][6 ];
        double clcdmin      = context.aerodata[n][7 ];
        double dcddcl2      = context.aerodata[n][8 ];
        double cmcon        = context.aerodata[n][9 ];
        double reref        = context.aerodata[n][10];
        double rexp         = context.aerodata[n][11];
        double mcrit        = context.aerodata[n][12];
        double xisect1      = context.xiaero[n];
        // Get data for inner bounding aero section
        clcdcm(context,
               alf, w, rey,
               clift, cl_alf, cl_w, stallf,
               cdrag, cd_alf, cd_w, cd_rey,
               cmom, cm_al, cm_w,
               a0, clmax, clmin, dclda, dclda_stall, dcl_stall,
               cdmin, clcdmin, dcddcl2, cmcon, mcrit, reref, rexp);

        // Check for another bounding section, if not we are done,
        // if we have another section linearly interpolate data to station IS
        if (n < context.naero-1) {
            double xisect2 = context.xiaero[n+1];
            double frac = (context.xi[is] - xisect1) / (xisect2 - xisect1);

            a0          = context.aerodata[n+1][0 ];
            double clmax2       = context.aerodata[n+1][1 ];
            double clmin2       = context.aerodata[n+1][2 ];
            dclda               = context.aerodata[n+1][3 ];
            dclda_stall         = context.aerodata[n+1][4 ];
            double dcl_stall2   = context.aerodata[n+1][5 ];
            cdmin               = context.aerodata[n+1][6 ];
            clcdmin             = context.aerodata[n+1][7 ];
            dcddcl2             = context.aerodata[n+1][8 ];
            cmcon               = context.aerodata[n+1][9 ];
            reref               = context.aerodata[n+1][10];
            rexp                = context.aerodata[n+1][11];
            mcrit               = context.aerodata[n+1][12];

            // Get data for outer bounding aero section
            double clift2, cl_alf2, cl_w2,
                    cdrag2, cd_alf2, cd_w2, cd_rey2,
                    cmom2, cm_al2, cm_w2;
            bool stallf2;
            clcdcm(context,
                   alf, w, rey,
                   clift2, cl_alf2, cl_w2, stallf2,
                   cdrag2, cd_alf2, cd_w2, cd_rey2,
                   cmom2, cm_al2, cm_w2,
                   a0, clmax2, clmin2, dclda, dclda_stall, dcl_stall2,
                   cdmin, clcdmin, dcddcl2, cmcon, mcrit, reref, rexp);

            // Interpolate aero data to blade station
            stallf = (stallf || stallf2);
            clift  = (1.0-frac)*clift  + frac*clift2;
            cl_alf = (1.0-frac)*cl_alf + frac*cl_alf2;
            cl_w   = (1.0-frac)*cl_w   + frac*cl_w2;
            clmax  = (1.0-frac)*clmax  + frac*clmax2;
            clmin  = (1.0-frac)*clmin  + frac*clmin2;
            dcl_stall = (1.0-frac)*dcl_stall + frac*dcl_stall2;

            cmom   = (1.0-frac)*cmom   + frac*cmom2;
            cm_al  = (1.0-frac)*cm_al  + frac*cm_al2;
            cm_w   = (1.0-frac)*cm_w   + frac*cm_w2;

            cdrag  = (1.0-frac)*cdrag  + frac*cdrag2;
            cd_alf = (1.0-frac)*cd_alf + frac*cd_alf2;
            cd_w   = (1.0-frac)*cd_w   + frac*cd_w2;
            cd_rey = (1.0-frac)*cd_rey + frac*cd_rey2;
        }
    }

    /**
     * Get the angle of attack for a given lift coefficient at blade station with index is.
     *
     * @param context
     * @param is        blade station index
     * @param clift     lift coefficient
     * @param w         ?
     * @param alf       angle of attack, alpha
     * @param alf_cl    d(alpha)/d(C_l)
     * @param alf_w     d(alpha)/d(w)
     * @param stallf    true if stalled
     */
    void getalf(common::context &context,
                unsigned is, const double &clift, const double &w,
                double &alf, double &alf_cl, double &alf_w, bool &stallf) {
        const int niter = 10;
        const double eps = 1.0e-5;

        stallf = false;

        // HHY had to set A0 to first aero section as A0 is now section property
        double a0 = context.aerodata[0][0];
        double rey = 0;

        alf = a0;
        double cltemp, cl_alf, cl_w, clmax, clmin, dcl_stall, cdrag, cd_alf, cd_w, cd_rey, cmom, cm_al, cm_w, dalf;
        for (int iter = 0; iter < niter; iter++) {
            getclcdcm(context,
                      is, alf, w, rey,
                      cltemp, cl_alf, cl_w,
                      clmax, clmin, dcl_stall, stallf,
                      cdrag, cd_alf, cd_w, cd_rey,
                      cmom, cm_al, cm_w);

            dalf = -(cltemp - clift) / cl_alf;
            alf += dalf;
            alf_cl =  1.0 / cl_alf;
            alf_w = -cl_w / cl_alf;
            if (abs(dalf) < eps) return;
        }
    }

    /**
     * Calculate C_l(alpha), C_d(alpha), and C_m(alpha).
     *
     * @see[xaero::putaero(), xaero::getclcdcm()] for description of input parameters.
     * 
     * #CL(alpha) function
     * Note that in addition to setting CLIFT and its derivatives
     * CLMAX and CLMIN (+ and - stall CL's) are set in this routine
     * In the compressible range the stall CL is reduced by a factor
     * proportional to Mcrit-Mach.  Stall limiting for compressible
     * cases begins when the compressible drag added CDC > CDMstall
     *
     * #CD(alpha) function - presently CD is assumed to be a sum
     * of profile drag + stall drag + compressibility drag
     * In the linear lift range drag is CD0 + quadratic function of CL-CLDMIN
     * In + or - stall an additional drag is added that is proportional
     * to the extent of lift reduction from the linear lift value.
     * Compressible drag is based on adding drag proportional to
     * (Mach-Mcrit_eff)^MEXP
     *
     * #CM(alpha) function - presently CM is assumed constant,
     * varying only with Mach by Prandtl-Glauert scaling
     */
    void clcdcm(common::context &context,
                const double &alf, const double &w, const double &rey,
                double &clift, double &cl_alf, double &cl_w, bool &stallf,
                double &cdrag, double &cd_alf, double &cd_w, double &cd_rey,
                double &cmom, double &cm_al, double &cm_w,
                const double &a0, double &clmax, double &clmin,
                const double &dclda, const double &dclda_stall, const double &dcl_stall,
                const double &cdmin, const double &cldmin, const double &dcddcl2,
                const double &cmcon, const double &mcrit, const double &reref, const double &rexp) {

        // Factors for compressibility drag model, HHY 10/23/00
        // Mcrit is set by user
        // Effective Mcrit is Mcrit_eff = Mcrit - CLMFACTOR*(CL-CLDmin) - DMDD
        // DMDD is the delta Mach to get CD=CDMDD (usually 0.0020)
        // Compressible drag is CDC = CDMFACTOR*(Mach-Mcrit_eff)^MEXP
        // CDMstall is the drag at which compressible stall begins
        const double cdmfactor = 10.0;
        const double clmfactor =  0.25;
        const double mexp      =  3.0;
        const double cdmstall  =  0.1000;

        // Prandtl-Glauert compressibility factor
        double msq   =   w * w * pow(context.vel, 2) / pow(context.vso, 2);
        double msq_w = 2.0 * w * pow(context.vel, 2) / pow(context.vso, 2);
        if (msq >= 1.0) {
            cout << "CLFUNC: Local Mach number limited to 0.99, was " << msq << endl;
            msq = 0.99;
            msq_w = 0;
        }
        double pg   = 1.0 / sqrt(1.0 - msq);
        double pg_w = 0.5 * msq_w * pow(pg, 3);

        // Mach number and dependence on velocity
        double mach     = sqrt(msq);
        double mach_w   = 0.0;
        if (mach != 0.0) mach_w = 0.5 * msq_w / mach;

        // =======================================================
        // Generate CL from dCL/dAlpha and Prandtl-Glauert scaling
        double cla      = dclda * pg   * (alf - a0);
        double cla_alf  = dclda * pg;
        double cla_w    = dclda * pg_w * (alf - a0);

        // Effective CLmax is limited by Mach effects
        // reduces CLmax to match the CL of onset of serious compressible drag
        double dmstall  = pow(cdmstall / cdmfactor, 1.0 / mexp);
        double clmaxm   = max(0.0, (mcrit + dmstall - mach) / clmfactor) + cldmin;
        clmax           = min(clmax, clmaxm);
        double clminm   = min(0.0, -(mcrit + dmstall - mach) / clmfactor) + cldmin;
        clmin           = max(clmin, clminm);

        // CL limiter function (turns on after +-stall
        double ecmax = exp( min(200.0, (cla - clmax) / dcl_stall) );
        double ecmin = exp( min(200.0, (clmin - cla) / dcl_stall) );
        double cllim = dcl_stall * log( (1.0 + ecmax) / (1.0 + ecmin) );
        double cllim_cla = ecmax / (1.0 + ecmax) + ecmin / (1.0 + ecmin);

        // Subtract off a (nearly unity) fraction of the limited CL function
        // This sets the dCL/dAlpha in the stalled regions to 1-FSTALL of that
        // in the linear lift range
        double fstall = dclda_stall / dclda;
        clift   = cla       - (1.0 - fstall) * cllim;
        cl_alf  = cla_alf   - (1.0 - fstall) * cllim_cla * cla_alf;
        cl_w    = cla_w     - (1.0 - fstall) * cllim_cla * cla_w;

        stallf = false;
        if (clift > clmax) stallf = true;
        if (clift < clmin) stallf = true;
        // =======================================================

        // =========================================
        // CM from CMCON and Prandtl-Glauert scaling
        cmom    = pg * cmcon;
        cm_al   = 0.0;
        cm_w    = pg * cmcon;
        // =========================================

        // =========================================================
        // CD from profile drag, stall drag and compressibility drag

        // Reynolds number scaling factor
        double rcorr, rcorr_rey;
        if (rey <= 0) {
            rcorr       = 1.0;
            rcorr_rey   = 0.0;
        } else {
            rcorr       = pow(rey / reref, rexp);
            rcorr_rey   = rexp / rey;
        }

        // In the basic linear lift range drag is a function of lift
        // CD = CD0 (constant) + quadratic with CL)
        cdrag   = (cdmin + dcddcl2 * pow(clift - cldmin, 2)      ) * rcorr;
        cd_alf  = (     2.0 * dcddcl2 * (clift - cldmin) * cl_alf) * rcorr;
        cd_w    = (     2.0 * dcddcl2 * (clift - cldmin) * cl_w  ) * rcorr;
        cd_rey  = cdrag * rcorr_rey;

        // Post-stall drag added
        fstall = dclda_stall / dclda;
        double dcdx     = (1.0 - fstall) * cllim / (pg * dclda);
        double dcd      = 2.0 * pow(dcdx, 2);
        double dcd_alf  = 4.0 * dcdx *
                         (1.0 - fstall) * cllim_cla * cla_alf / (pg * dclda);
        double dcd_w    = 4.0 * dcdx *
                       ( (1.0 - fstall) * cllim_cla * cla_w   / (pg * dclda) - dcd / pg * pg_w );

        // Compressibility drag (accounts for drag rise above Mcrit with CL effects
        // CDC is a function of a scaling factor*(M-Mcrit(CL))**MEXP
        // DMDD is the Mach difference corresponding to CD rise of CDMDD at MCRIT
        double critmach     = mcrit - clmfactor * abs(clift);
        double critmach_alf = -clmfactor * abs(cl_alf);
        double critmach_w   = -clmfactor * abs(cl_w);
        double cdc, cdc_alf, cdc_w;
        if (mach < critmach) {
            cdc     = 0.0;
            cdc_alf = 0.0;
            cdc_w   = 0.0;
        } else {
            cdc     = cdmfactor * pow(mach - critmach, mexp);
            cdc_alf =                            - mexp * critmach_alf * cdc / critmach;
            cdc_w   = mexp * mach_w * cdc / mach - mexp * critmach_w   * cdc / critmach;
        }

        double fac   = 1.0;
        double fac_w = 0.0;
        // Although test data does not show profile drag increases due to Mach #
        // you could use something like this to add increase drag by Prandtl-Glauert
        // (or any function you choose)
        // fac      = pg;
        // fac_w    = pg_w;

        // Total drag terms
        cdrag   = fac * cdrag                + dcd     + cdc;
        cd_alf  = fac * cd_alf               + dcd_alf + cdc_alf;
        cd_w    = fac * cd_w + fac_w * cdrag + dcd_w   + cdc_w;
        cd_rey  = fac * cd_rey;
        // =========================================================
    }
}