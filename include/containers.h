//
// Created by daniel.devries on 11/27/2018.
//

#ifndef XROTOR_NOGRAPHICS_C_ARRAYS_H
#define XROTOR_NOGRAPHICS_C_ARRAYS_H

class Matrix {
public:
    Matrix(unsigned long m, unsigned long n);
    Matrix(unsigned long m, unsigned long n, double val);
    Matrix(const Matrix& mat);
    ~Matrix();

    double& operator() (unsigned long i, unsigned long j);
    double operator() (unsigned long i, unsigned long j) const;

    Matrix& zeros();

private:
    unsigned long m_, n_;
    double* data_;
};

class Cube {
public:
    Cube(unsigned long m, unsigned long n, unsigned long o);
    Cube(unsigned long m, unsigned long n, unsigned long o, double val);
    Cube(const Cube& cube);
    ~Cube();

    double& operator() (unsigned long i, unsigned long j, unsigned long k);
    double operator() (unsigned long i, unsigned long j, unsigned long k) const;

    Cube& zeros();

private:
    unsigned long m_, n_, o_;
    double* data_;
};

#endif //XROTOR_NOGRAPHICS_C_ARRAYS_H
