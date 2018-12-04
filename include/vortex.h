//
// Created by Daniël de Vries on 26/11/2018.
//

#ifndef XROTOR_NOGRAPHICS_CPP_VORTEX_H
#define XROTOR_NOGRAPHICS_CPP_VORTEX_H

namespace vortex {
    template <int imax>
    void vrtxc0(int ii, int nblds, bool lduct, double rake,
                const double xi[imax], const double xv[imax], const double gam[imax], double adw,
                double vind_gam[3][imax][imax], double vind_adw[3][imax]);

    void vorsegvel(const double a[3], const double b[3], double uvw[3], double uvw_a[3][3], double uvw_b[3][3]);
}

#endif //XROTOR_NOGRAPHICS_CPP_VORTEX_H
