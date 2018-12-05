//
// Created by daniel.devries on 11/29/2018.
//

#ifndef XROTOR_NOGRAPHICS_C_XIO_H
#define XROTOR_NOGRAPHICS_C_XIO_H

#include <fstream>      // std::ifstream
#include <string>       // std::string
using namespace std;

#include <common.h>

namespace xio {
    void LOAD(common::context &ctxt, string fname1);
    void RDLINE(ifstream &ifs, string &line);
    bool READ(ifstream &ifs, string &line, const char *fmt, int n, ...);
}

#endif //XROTOR_NOGRAPHICS_C_XIO_H
