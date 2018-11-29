//
// Created by Daniël de Vries on 26/11/2018.
//

#ifndef XROTOR_NOGRAPHICS_CPP_COMMON_H
#define XROTOR_NOGRAPHICS_CPP_COMMON_H

#include <vector>
using std::vector;
typedef std::vector<double> vec;
typedef std::vector<unsigned> uvec;

namespace common {
    const int ix = 100;                 /**< max number of radial prop stations */
    const int ipx = ix + 1;
    const int nparx = 12;               /**< max number of stored cases */
    const int icasx = 100;              /**< number of case parameters stored */
    const int iwx = 200;                /**< dimension of work arrays */

    const int nax = 20;                 /**< max number of aerodynamic sections defined */
    const int ndx = 14;                 /**< number of aerodynamic parameters defined for each section */

    const int iq = ix + 5;
    const int jx = (iq * 3) / 2 + 1;

    const double pi = 3.14159265358979323846264338327950288;

    struct context {
        double mcrit;
        double rho, rmu, vso, vel, rad, gee, alt;

        unsigned long ii;
        unsigned naero;

        uvec iaero;
        vec xiaero, xi;
        vector<vector<double>> aerodata;

    };
}






#endif //XROTOR_NOGRAPHICS_CPP_COMMON_H
