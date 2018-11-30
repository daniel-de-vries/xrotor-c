//
// Created by daniel.devries on 11/30/2018.
//

#ifndef XROTOR_NOGRAPHICS_C_XROTOR_H
#define XROTOR_NOGRAPHICS_C_XROTOR_H

#include <fstream>
#include "common.h"

using namespace std;

namespace xrotor {
    void xrotor();
    void init(common::context& context);
    void setdef(common::context& context);
    void atmo(double alspec, double& vsoalt, double& rhoalt, double& rmualt);
    void flosho(double vso, double rho, double rmu);
    void reinit(common::context& context);
    void setx(common::context& context);
    void opfile(ofstream& ofs, string& fname);
    void output(common::context &context, ostream &ofs);
    void uvadd(common::context& context, double xiw, double& wa, double& wt);
}

#endif //XROTOR_NOGRAPHICS_C_XROTOR_H
