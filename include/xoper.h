//
// Created by daniel.devries on 11/30/2018.
//

#ifndef XROTOR_NOGRAPHICS_C_XOPER_H
#define XROTOR_NOGRAPHICS_C_XOPER_H

#include <vector>
using namespace std;

#include <common.h>

namespace xoper {
    void OPER(common::context &context);

    namespace commands {
        void FORM(common::context& context);
        void TERS(common::context& context);
        void DISP(common::context& context);
        void NAME(common::context& context);
        void WRIT(common::context& context, string fName);
        void DUCT(common::context& context);
        void VRAT(common::context& context);
        void ATMO(common::context& context);
        void VELO(common::context& context);
        void ANGL(common::context& context);
        void ADVA(common::context& context);
        void RPM (common::context& context, int ninput, vector<double> rinput);
        void THRU(common::context& context, int ninput, vector<double> rinput);
        void TORQ(common::context& context);
        void POWE(common::context& context);
        void ASEQ(common::context& context);
        void RSEQ(common::context& context);
        void BSEQ(common::context& context);
        void VSEQ(common::context& context);
        void CLRC(common::context& context);
        void ADDC(common::context& context);
        void CPUT(common::context& context);
        void CGET(common::context& context);
        void CASE(common::context& context);
        void LIST(common::context& context);
        void N   (common::context& context);
        void ITER(common::context& context, int ninput, vector<int> iinput);
        void INIT(common::context& context);
        void REIN(common::context& context);
    }

    void APER  (common::context& context, unsigned short ispec, unsigned short icon, bool linit);
    void APINIT(common::context& context);
    void APITER(common::context& context, unsigned short ispec, bool icon);

    void CSCALC(common::context& context,
                const int &i, double& utot, const double &wa, const double &wt,
                double & vt, double &vt_adw,
                double & va, double &va_adw,
                double & vd, double &vd_adw,
                double & ci, double &ci_adv, double &ci_vt,
                double & si,                                double &si_va,
                double &  w, double & w_adv, double & w_vt, double & w_va,
                double &phi, double & p_adv, double & p_vt, double & p_va);

    void XWINIT(common::context& context);
    void SETXW (common::context& context);
    void TPQ   (common::context& context, unsigned short itype);
    void VCALC (common::context& context);

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
