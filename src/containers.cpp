//
// Created by daniel.devries on 11/28/2018.
//

#include <algorithm>    // std::copy
#include <stdexcept>    // std::out_of_range
#include "containers.h"

// Matrix DEFINITION START
inline
Matrix::Matrix(unsigned m, unsigned n)
        : m_ (m)
        , n_ (n)
{
    if (m == 0 || n == 0)
        throw std::out_of_range("m and n should be > 0");
    data_ = new double[m * n];
}

inline
Matrix::Matrix(const Matrix &mat) {
    m_ = mat.m_;
    n_ = mat.n_;
    unsigned N = m_ * n_;
    data_ = new double[N];
    std::copy(mat.data_, mat.data_ + N, data_);
}

inline
Matrix::~Matrix() {
    delete[] data_;
}

inline
double& Matrix::operator()(unsigned i, unsigned j) {
    if (i >= m_ || j >= n_)
        throw std::out_of_range("Indices out of bounds.");
    return data_[i*n_ + j];
}

inline
double Matrix::operator()(unsigned i, unsigned j) const {
    if (i >= m_ || j >= n_)
        throw std::out_of_range("Indices out of bounds.");
    return data_[i*n_ + j];
}
// Matrix DEFINITION END

// Cube DEFINITION START
inline
Cube::Cube(unsigned m, unsigned n, unsigned o)
        : m_ (m)
        , n_ (n)
        , o_ (o)
{
    if (m == 0 || n == 0 || o == 0)
        throw std::out_of_range("m, n, and o should be > 0");
    data_ = new double[m * n * o];
}

inline
Cube::Cube(const Cube &cube) {
    m_ = cube.m_;
    n_ = cube.n_;
    o_ = cube.o_;
    unsigned N = m_ * n_ * o_;
    data_ = new double[N];
    std::copy(cube.data_, cube.data_ + N, data_);
}

inline
Cube::~Cube() {
    delete[] data_;
}

inline
double& Cube::operator()(unsigned i, unsigned j, unsigned k) {
    if (i >= m_ || j >= n_ || k >= o_)
        throw std::out_of_range("Indices out of bounds.");
    return data_[i*n_*o_ + j*o_ + k];
}

inline
double Cube::operator()(unsigned i, unsigned j, unsigned k) const {
    if (i >= m_ || j >= n_ || k >= o_)
        throw std::out_of_range("Indices out of bounds.");
    return data_[i*n_*o_ + j*o_ + k];
}
// Cube DEFINITION END
