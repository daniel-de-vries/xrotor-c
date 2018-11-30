//
// Created by Daniël de Vries on 26/11/2018.
//

#ifndef XROTOR_NOGRAPHICS_CPP_XAERO_H
#define XROTOR_NOGRAPHICS_CPP_XAERO_H

#include "common.h"

namespace xaero {
    void setiaero(common::context &context);

    void putaero(common::context &context,
                 int n, double xisect, double a0, double clmax, double clmin,
                 double dclda, double dclda_stall, double dcl_stall,
                 double cdmin, double clcdmin, double dcddcl2,
                 double cmcon, double mcrit, double reref, double rexp);

    void getclcdcm(common::context &context,
                   int is, double alf, double w, double rey,
                   double clfit, double cl_alf, double cl_w,
                   double clmax, double clmin, double dcl_stall, double stallf,
                   double cdrag, double cd_alf, double cd_w, double cd_rey,
                   double cmom, double cm_al, double cm_w);

    void getalf(common::context &context,
                int is, const double &clift, const double &w,
                double &alf, double &alf_cl, double &alf_w, bool &stallf);

    void clcdcm(common::context &context,
                const double &alf, const double &w, const double &rey,
                double &clift, double &cl_alf, double &cl_w, bool &stallf,
                double &cdrag, double &cd_alf, double &cd_w, double &cd_rey,
                double &cmom, double &cm_al, double &cm_w,
                const double &a0, double &clmax, double &clmin,
                const double &dclda, const double &dclda_stall, const double &dcl_stall,
                const double &cdmin, const double &cldmin, const double &dcddcl2,
                const double &cmcon, const double &mcrit, const double &reref, const double &rexp);
}

#endif //XROTOR_NOGRAPHICS_CPP_XAERO_H
