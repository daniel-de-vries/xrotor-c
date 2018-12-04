//
// Created by daniel.devries on 11/30/2018.
//

#ifndef XROTOR_NOGRAPHICS_C_XUTILS_H
#define XROTOR_NOGRAPHICS_C_XUTILS_H

namespace xutils {
    template <int nsiz, int nrhs_full>
    void GAUSS(int nn, double z[nsiz][nsiz], double r[nsiz][nrhs_full], int nrhs);
}

#endif //XROTOR_NOGRAPHICS_C_XUTILS_H
