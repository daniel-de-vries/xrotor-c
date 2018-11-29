//
// Created by Daniël de Vries on 26/11/2018.
//

#ifndef XROTOR_NOGRAPHICS_CPP_SPLINE_H
#define XROTOR_NOGRAPHICS_CPP_SPLINE_H

#include <vector>
typedef std::vector<double> vec;

namespace spline {

    void spline(const vec &x, vec &xs, const vec &s);
    void splind(const vec &x, vec &xs, const vec &s, double xs1, double xs2);
    void splina(const vec &x, vec &xs, const vec &s);
    void trisol(vec &a, const vec &b, vec &c, vec &d);

    double seval(double ss, const vec &x, const vec &xs, const vec &s);
    double deval(double ss, const vec &x, const vec &xs, const vec &s);

    void segspl(const vec &x, vec &xs, const vec &s);
}

#endif // XROTOR_NOGRAPHICS_CPP_SPLINE_H
