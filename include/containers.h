//
// Created by daniel.devries on 11/27/2018.
//

#ifndef XROTOR_NOGRAPHICS_C_ARRAYS_H
#define XROTOR_NOGRAPHICS_C_ARRAYS_H

class Matrix {
public:
    Matrix(int m, int n);
    Matrix(int m, int n, double val);
    Matrix(const Matrix& mat);
    ~Matrix();

    double& operator() (int i, int j);
    double operator() (int i, int j) const;

    Matrix& zeros();

private:
     int m_, n_;
    double* data_;
};

class Cube {
public:
    Cube( int m,  int n,  int o);
    Cube( int m,  int n,  int o, double val);
    Cube(const Cube& cube);
    ~Cube();

    double& operator() ( int i,  int j,  int k);
    double operator() ( int i,  int j,  int k) const;

    Cube& zeros();

private:
     int m_, n_, o_;
    double* data_;
};

#endif //XROTOR_NOGRAPHICS_C_ARRAYS_H
