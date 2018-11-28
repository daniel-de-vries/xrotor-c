//
// Created by Daniël de Vries on 26/11/2018.
//
#include <algorithm>    // std::copy
#include <stdexcept>    // std::runtime_error
#include "spline.h"

namespace spline {

    /**
     * Calculate the spline coefficients for x(s).
     *
     * Zero 2nd derivative end conditions are used.
     * To evaluate the spline at some value of s,
     * use spline::seval() and/or spline::deval().
     *
     * @param s     independent variable array
     * @param x     dependent variable array
     * @param xs    d(x)/d(s) array (calculated)
     */
    void spline(const vec &x, vec &xs, const vec &s) {
        return spline::splind(x, xs, s, 999.f, 999.f);
    }

    /**
     * Calculate the spline coefficients for x(s).
     *
     * Specified 1st derivative and/or zero 2nd
     * or 3rd derivative end conditions can be used.
     * To evaludate the spline at some value of s,
     * use spline::seval() and/or spline::deval().
     *
     * @param s         independent variable array
     * @param x         dependent variable array
     * @param xs        d(x)/d(s) array (calculated)
     * @param xs1, xs2  endpoint derivatives
     *                  if ==  999.0, use zero 2nd derivative
     *                  if == -999.0, use zero 2rd derivative
     */
    void splind(const vec &x, vec &xs, const vec &s, double xs1, double xs2) {
        unsigned long n = x.size();
        if (n == 1) {
            xs[0] = .0f;
            return;
        }

        vec a(n), b(n), c(n);

        for (unsigned long i = 1; i < n -1; i++) {
            double dsm = s[i] - s[i-1];
            double dsp = s[i+1] - s[i];
            b[i] = dsp;
            a[i] = 2 * (dsm + dsp);
            c[i] = dsm;

            xs[i] = 3 * ((x[i+1] - x[i]) * dsm / dsp + (x[i] - x[i-1]) * dsp / dsm);
        }

        // set left end condition
        if (xs1 == 999.f) {
            // set zero 2nd derivative
            a[0] = 2.f;
            c[0] = 1.f;
            xs[0] = 3 * (x[1] - x[0]) / (s[1] - s[0]);
        } else if (xs1 == -999.f) {
            // set zero 3rd derivative
            a[0] = 1.f;
            c[0] = 1.f;
            xs[0] = 2 * (x[1] - x[0]) / (s[1] - s[0]);
        } else {
            // set specified 1st derivative
            a[0] = 1.f;
            c[0] = 0.f;
            xs[0] = xs1;
        }

        // set right conditions
        if (xs2 == 999.f) {
            // set zero 2nd derivative
            a[n-1] = 2.f;
            b[n-1] = 1.f;
            xs[n-1] = 3 * (x[n-1] - x[n-2]) / (s[n-1] - s[n-2]);
        } else if (xs2 == -999.f) {
            // set zero 3rd derivative
            a[n-1] = 1.f;
            b[n-1] = 1.f;
            xs[n-1] = 2 * (x[n-1] - x[n-2]) / (s[n-1] - s[n-2]);
        } else {
            // set specified 1st derivative
            a[n-1] = 1.f;
            b[n-1] = 0.f;
            xs[n-1] = xs2;
        }

        // solve for derivative array xs
        return spline::trisol(a, b, c, xs);
    }

    /**
     * Calculate the spline coefficients for x(s).
     *
     * A simple averaging of adjecent segment slopes
     * is used to achieve non-oscillatory curve.
     * End conditions are set by end segment slope.
     * To evaluate the spline at some value of s,
     * use spline::seval() and/or spline::deval().
     *
     * @param s     independent variable array
     * @param x     dependent variable array
     * @param xs    d(x)/d(s) array
     *
     * @note All three arrays should be conformant.
     *
     */
    void splina(const vec &x, vec &xs, const vec &s) {
        unsigned long n = x.size();
        if (n == 1) {
            xs[0] = 0.f;
            return;
        }

        bool lend = true;
        double xs1 = 0;
        double xs2 = 0;
        for (unsigned long long i = 0; i < n-2; i++) {
            double ds = s[i+1] - s[i];

            if (ds == 0.f) {
                xs[i] = xs1;
            } else {
                double dx = x[i+1] - x[i];
                xs2 = dx / ds;
                if (lend) {
                    xs[i] = xs2;
                    lend = false;
                } else {
                    xs[i] = .5f * (xs1 + xs2);
                }
            }

            xs1 = xs2;
        }
        xs[n-1] = xs1;
    }

    /**
     * Solve long, tri-diagonal system
     *
     *         a c          d
     *         b a c        d
     *           b a .      .
     *             . . c    .
     *               b a    d
     *
     * The right hand side, d, is replaced by
     * the solution. a and c are destroyed.
     *
     * @param a, b, c, d    arrays
     */
    void trisol(vec &a, const vec &b, vec &c, vec &d) {
        unsigned long kk = a.size();
        for (unsigned long k = 1; k < kk-1; k++) {
            c[k-1] /= a[k-1];
            d[k-1] /= a[k-1];
            a[k] -= b[k] * c[k-1];
            d[k] -= b[k] * d[k-1];
        }

        d[kk-1] /= a[kk-1];

        for (unsigned long k = kk-2; k >= 0; k--) {
            d[k] -= c[k] * d[k+1];
        }
    }

    /**
     * Helper function for spline::seval() and
     * spline::deval() to avoid code duplication.
     *
     * @see[spline::seval, spline::deval].
     */
    void _eval_helper(double ss, const vec &x, const vec &xs, const vec &s,
                      unsigned long &i, double &ds, double &t, double &cx1, double &cx2) {
        unsigned long ilow = 0;
        i = x.size()-1;
        while (i - ilow > 1) {
            unsigned long imid = (i + ilow) / 2;
            if (ss < s[imid]) {
                i = imid;
            } else {
                ilow = imid;
            }
        }

        ds = s[i] - s[i-1];
        t = (ss - s[i-1]) / ds;
        cx1 = ds * xs[i-1] - x[i] + x[i-1];
        cx2 = ds * xs[i]   - x[i] + x[i-1];
    }

    /**
     * Calculate x(ss).
     *
     * xs array must have been calculated by spline::spline().
     *
     * @param ss    point at which to evaluate the spline
     * @param s     independent variable array
     * @param x     dependent variable array
     * @param xs    d(x)/d(s) array
     * @return      x(ss)
     */
    double seval(double ss, const vec &x, const vec &xs, const vec &s) {
        unsigned long n = x.size();
        if (n == 1)
            return xs[0];

        unsigned long i;
        double ds, t, cx1, cx2;
        spline::_eval_helper(ss, x, xs, s, i, ds, t, cx1, cx2);
        return t * x[i] + (1 - t) * x[i-1] + (t - t*t) * ((1 - t) * cx1 - t * cx2);
    }

    /**
     * Calculate [d(x)d(s)](ss).
     *
     * xs array must have been calculated by spline::spline().
     *
     * @param ss    point at which to evaluate the spline
     * @param s     independent variable array
     * @param x     dependent variable array
     * @param xs    d(x)/d(s) array
     * @return      [d(x)/d(s)](ss)
     */
    double deval(double ss, const vec &x, const vec &xs, const vec &s) {
        unsigned long n = x.size();
        if (n == 1)
            return xs[0];

        unsigned long i;
        double ds, t, cx1, cx2;
        spline::_eval_helper(ss, x, xs, s, i, ds, t, cx1, cx2);
        return (x[i] - x[i-1] + (1 - 4*t + 3*t*t) * cx1 + t * (3*t - 2) * cx2) / ds;
    }

    /**
     * Helper function for spline::segspl() to avoid code duplication.
     *
     * @see[spline::segspl].
     */
    void _segspl_helper(const vec &x, vec &xs, const vec &s, unsigned long iseg, unsigned long iseg0) {
        unsigned long nseg = iseg - iseg0 + 1;
        vec _x(nseg), _xs(nseg), _s(nseg);
        std::copy(x.begin() + iseg0, x.end(), _x.begin());
        std::copy(xs.begin() + iseg0, xs.end(), _xs.begin());
        std::copy(s.begin() + iseg0, s.end(), _s.begin());
        spline::splind(_x, _xs, _s, -999.f, -999.f);
        std::copy(_xs.begin(), _xs.end(), xs.begin() + iseg0);
    }

    /**
     * Spline x(s) array just like spline::spline(),
     * but allow derivative discontinuities
     * at segment joins. Segment joints are
     * defined by identical successive s values.
     *
     * @param s     independent variable array
     * @param x     dependent variable array
     * @param xs    d(x)/d(s) array (calculated)
     */
    void segspl(const vec &x, vec &xs, const vec &s) {
        unsigned long n = x.size();
        if (n == 1) {
            xs[0] = 0.f;
            return;
        }

        if (s[0]   == s[1]  ) throw std::runtime_error("segspl:  First input point duplicated");
        if (s[n-1] == s[n-2]) throw std::runtime_error("segspl:  Last  input point duplicated");

        unsigned long iseg0 = 0;
        for (unsigned long iseg = 1; iseg < n-3; iseg++) {
            if (s[iseg] == s[iseg+1]) {
                _segspl_helper(x, xs, s, iseg, iseg0);
                iseg0 = iseg + 1;
            }
        }

        _segspl_helper(x, xs, s, n, iseg0);
    }
}

