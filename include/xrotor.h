//
// Created by daniel.devries on 11/30/2018.
//

#ifndef XROTOR_NOGRAPHICS_C_XROTOR_H
#define XROTOR_NOGRAPHICS_C_XROTOR_H

#include <fstream>

#include <common.h>

using namespace std;

namespace xrotor {
    void XROTOR();
    void INIT(common::context &context);
    void SETDEF(common::context &context);
    void ATMO(double alspec, double &vsoalt, double &rhoalt, double &rmualt);
    void FLOSHO(ostream &os, double vso, double rho, double rmu);
    void REINIT(common::context &context);
    void SETX(common::context &context);
    void OPFILE(ofstream &ofs, string &fname);
    void OUTPUT(common::context &context, ostream &ofs);
    void UVADD(common::context &context, double xiw, double &wa, double &wt);
}

#endif //XROTOR_NOGRAPHICS_C_XROTOR_H
