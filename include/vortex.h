//
// Created by Daniël de Vries on 26/11/2018.
//

#ifndef XROTOR_NOGRAPHICS_CPP_VORTEX_H
#define XROTOR_NOGRAPHICS_CPP_VORTEX_H

#include "containers.h" // Matrix, Cube
#include <vector>       // std::vector
typedef std::vector<double> vec;

namespace vortex {
    void vrtxc0(int ii, int nblds, bool lduct, double rake,
                const vec &xi, const vec &xv, const vec &gam, double adw,
                Cube &vind_gam, Matrix &vind_adw);

    void vorsegvel(const double a[3], const double b[3], double uvw[3], Matrix &uvw_a, Matrix &uvw_b);
}

#endif //XROTOR_NOGRAPHICS_CPP_VORTEX_H
