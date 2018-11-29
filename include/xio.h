//
// Created by daniel.devries on 11/29/2018.
//

#ifndef XROTOR_NOGRAPHICS_C_XIO_H
#define XROTOR_NOGRAPHICS_C_XIO_H

#include <fstream>      // std::ifstream
#include <string>       // std::string
#include "common.h"
using namespace std;

namespace xio {
    void load(common::context &context, string fname1);
    void rdline(ifstream &ifs, string &line);
    bool read(ifstream &ifs, string &line, const char *fmt, int n, ...);
}

#endif //XROTOR_NOGRAPHICS_C_XIO_H
