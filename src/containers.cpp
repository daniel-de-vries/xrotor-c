//
// Created by daniel.devries on 11/28/2018.
//

#include <algorithm>    // std::copy
#include <stdexcept>    // std::out_of_range
#include "containers.h"

// Matrix DEFINITION START
Matrix::Matrix(unsigned long m, unsigned long n)
        : m_ (m)
        , n_ (n)
{
    if (m == 0 || n == 0)
        throw std::out_of_range("m and n should be > 0");
    data_ = new double[m * n];
}

Matrix::Matrix(unsigned long m, unsigned long n, double val) : Matrix::Matrix(m, n) {
    for (unsigned long i = 0; i < m*n; i++) {
        data_[i] = val;
    }
}

Matrix::Matrix(const Matrix &mat) {
    m_ = mat.m_;
    n_ = mat.n_;
    unsigned long N = m_ * n_;
    data_ = new double[N];
    std::copy(mat.data_, mat.data_ + N, data_);
}

Matrix::~Matrix() {
    delete[] data_;
}

double& Matrix::operator()(unsigned long i, unsigned long j) {
    if (i >= m_ || j >= n_)
        throw std::out_of_range("Indices out of bounds.");
    return data_[i*n_ + j];
}

double Matrix::operator()(unsigned long i, unsigned long j) const {
    if (i >= m_ || j >= n_)
        throw std::out_of_range("Indices out of bounds.");
    return data_[i*n_ + j];
}

Matrix& Matrix::zeros() {
    for (unsigned long i = 0; i < m_*n_; i++) {
        data_[i] = 0.;
    }
    return *this;
}
// Matrix DEFINITION END

// Cube DEFINITION START
Cube::Cube(unsigned long m, unsigned long n, unsigned long o)
        : m_ (m)
        , n_ (n)
        , o_ (o)
{
    if (m == 0 || n == 0 || o == 0)
        throw std::out_of_range("m, n, and o should be > 0");
    data_ = new double[m * n * o];
}

Cube::Cube(unsigned long m, unsigned long n, unsigned long o, double val) : Cube::Cube(m, n, o) {
    for (unsigned long i = 0; i < m*n*o; i++) {
        data_[i] = val;
    }
}

Cube::Cube(const Cube &cube) {
    m_ = cube.m_;
    n_ = cube.n_;
    o_ = cube.o_;
    unsigned long N = m_ * n_ * o_;
    data_ = new double[N];
    std::copy(cube.data_, cube.data_ + N, data_);
}

Cube::~Cube() {
    delete[] data_;
}

double& Cube::operator()(unsigned long i, unsigned long j, unsigned long k) {
    if (i >= m_ || j >= n_ || k >= o_)
        throw std::out_of_range("Indices out of bounds.");
    return data_[i*n_*o_ + j*o_ + k];
}

double Cube::operator()(unsigned long i, unsigned long j, unsigned long k) const {
    if (i >= m_ || j >= n_ || k >= o_)
        throw std::out_of_range("Indices out of bounds.");
    return data_[i*n_*o_ + j*o_ + k];
}

Cube& Cube::zeros() {
    for (unsigned long i = 0; i < m_*n_*o_; i++) {
        data_[i] = 0.;
    }
    return *this;
}
// Cube DEFINITION END
