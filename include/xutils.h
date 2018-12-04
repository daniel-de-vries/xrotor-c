//
// Created by daniel.devries on 11/30/2018.
//

#ifndef XROTOR_NOGRAPHICS_C_XUTILS_H
#define XROTOR_NOGRAPHICS_C_XUTILS_H

namespace xutils {
    template <int nsiz, int nrhs>
    void GAUSS(int nn, double z[nsiz][nsiz], double r[nsiz][nrhs]);
}

#endif //XROTOR_NOGRAPHICS_C_XUTILS_H
