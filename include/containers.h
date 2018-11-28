//
// Created by daniel.devries on 11/27/2018.
//

#ifndef XROTOR_NOGRAPHICS_C_ARRAYS_H
#define XROTOR_NOGRAPHICS_C_ARRAYS_H

#include <vector>

class Matrix {
public:
    Matrix(unsigned m, unsigned n);
    Matrix(const Matrix& mat);
    ~Matrix();

    double& operator() (unsigned i, unsigned j);
    double operator() (unsigned i, unsigned j) const;

private:
    unsigned m_, n_;
    double* data_;
};

class Cube {
public:
    Cube(unsigned m, unsigned n, unsigned o);
    Cube(const Cube& cube);
    ~Cube();

    double& operator() (unsigned i, unsigned j, unsigned k);
    double operator() (unsigned i, unsigned j, unsigned k) const;

private:
    unsigned m_, n_, o_;
    double* data_;
};

#endif //XROTOR_NOGRAPHICS_C_ARRAYS_H
