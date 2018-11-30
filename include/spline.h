//
// Created by Daniël de Vries on 26/11/2018.
//

#ifndef XROTOR_NOGRAPHICS_CPP_SPLINE_H
#define XROTOR_NOGRAPHICS_CPP_SPLINE_H

namespace spline {

    void spline(const double x[], double xs[], const double s[], int n);
    void splind(const double x[], double xs[], const double s[], int n, double xs1, double xs2);
    void splina(const double x[], double xs[], const double s[], int n);
    void trisol(double a[], const double b[], double c[], double d[], int kk);

    double seval(double ss, const double x[], double xs[], const double s[], int n);
    double deval(double ss, const double x[], double xs[], const double s[], int n);

    void segspl(const double x[], double xs[], const double s[], int n);
}

#endif // XROTOR_NOGRAPHICS_CPP_SPLINE_H
