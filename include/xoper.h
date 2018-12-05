//
// Created by daniel.devries on 11/30/2018.
//

#ifndef XROTOR_NOGRAPHICS_C_XOPER_H
#define XROTOR_NOGRAPHICS_C_XOPER_H

#include <vector>
using namespace std;

#include <common.h>

namespace xoper {
    void OPER(common::context &ctxt);

    namespace commands {
        void FORM(common::context &ctxt);
        void TERS(common::context &ctxt);
        void DISP(common::context &ctxt);
        void NAME(common::context &ctxt);
        void WRIT(common::context &ctxt, string fName);
        void DUCT(common::context &ctxt);
        void VRAT(common::context &ctxt);
        void ATMO(common::context &ctxt);
        void VELO(common::context &ctxt);
        void ANGL(common::context &ctxt);
        void ADVA(common::context &ctxt);
        void RPM (common::context &ctxt, int ninput, vector<double> rinput);
        void THRU(common::context &ctxt, int ninput, vector<double> rinput);
        void TORQ(common::context &ctxt);
        void POWE(common::context &ctxt);
        void ASEQ(common::context &ctxt);
        void RSEQ(common::context &ctxt);
        void BSEQ(common::context &ctxt);
        void VSEQ(common::context &ctxt);
        void CLRC(common::context &ctxt);
        void ADDC(common::context &ctxt);
        void CPUT(common::context &ctxt);
        void CGET(common::context &ctxt);
        void CASE(common::context &ctxt);
        void LIST(common::context &ctxt);
        void N   (common::context &ctxt);
        void ITER(common::context &ctxt, int ninput, vector<int> iinput);
        void INIT(common::context &ctxt);
        void REIN(common::context &ctxt);
    }

    void APER  (common::context &ctxt, unsigned short ispec, unsigned short icon, bool linit);
    void APINIT(common::context &ctxt);
    void APITER(common::context &ctxt, unsigned short ispec, bool icon);

    void CSCALC(common::context &ctxt,
                const int &i, double &utot, double &wa, double &wt,
                double & vt, double &vt_adw,
                double & va, double &va_adw,
                double & vd, double &vd_adw,
                double & ci, double &ci_adv, double &ci_vt,
                double & si,                                double &si_va,
                double &  w, double & w_adv, double & w_vt, double & w_va,
                double &phi, double & p_adv, double & p_vt, double & p_va);

    void XWINIT(common::context &ctxt);
    void SETXW (common::context &ctxt);
    void TPQ   (common::context &ctxt, unsigned short itype);
    void VCALC (common::context &ctxt);

    template <const int imax>
    void GRADMO(const int &ii, const int &nblds, const bool &lduct, const double &rake,
                const double xi[imax], const double xv[imax], const double gam[imax], const double &adw,
                double vind_gam[3][imax][imax], double vind_adw[3][imax]);

    template <const int imax>
    void HELICO(const int &ii, const int &nblds, const bool &lduct, const double &rake,
                const double xi[imax], const double xv[imax], const double gam[imax], const double &adw,
                double vind_gam[3][imax][imax], double vind_adw[3][imax]);

    void FILTER(double q[], const double &smlen, const int &n);

}

#endif //XROTOR_NOGRAPHICS_C_XOPER_H
