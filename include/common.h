//
// Created by Daniël de Vries on 26/11/2018.
//

#ifndef XROTOR_NOGRAPHICS_CPP_COMMON_H
#define XROTOR_NOGRAPHICS_CPP_COMMON_H

#include <armadillo>

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

    const double pi = arma::datum::pi;
}






#endif //XROTOR_NOGRAPHICS_CPP_COMMON_H
