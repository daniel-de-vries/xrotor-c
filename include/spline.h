//
// Created by Daniël de Vries on 26/11/2018.
//

#ifndef XROTOR_NOGRAPHICS_CPP_SPLINE_H
#define XROTOR_NOGRAPHICS_CPP_SPLINE_H

namespace spline {

    void SPLINE(const double *x, double *xs, const double *s, int n);
    void SPLIND(const double *x, double *xs, const double *s, int n, double xs1, double xs2);
    void SPLINA(const double *x, double *xs, const double *s, int n);
    void TRISOL(double *a, const double *b, double *c, double *d, int kk);

    double SEVAL(double ss, const double *x, double *xs, const double *s, int n);
    double DEVAL(double ss, const double *x, double *xs, const double *s, int n);

    void SEGSPL(const double *x, double *xs, const double *s, int n);
}

#endif // XROTOR_NOGRAPHICS_CPP_SPLINE_H
