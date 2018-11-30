//
// Created by daniel.devries on 11/28/2018.
//

#include <algorithm>    // std::copy
#include <stdexcept>    // std::out_of_range
#include "containers.h"

// Matrix DEFINITION START
Matrix::Matrix(int m,  int n)
        : m_ (m)
        , n_ (n)
{
    if (m == 0 || n == 0)
        throw std::out_of_range("m and n should be > 0");
    data_ = new double[m * n];
}

Matrix::Matrix(int m,  int n, double val) : Matrix::Matrix(m, n) {
    for (int i = 0; i < m*n; i++) {
        data_[i] = val;
    }
}

Matrix::Matrix(const Matrix &mat) {
    m_ = mat.m_;
    n_ = mat.n_;
     int N = m_ * n_;
    data_ = new double[N];
    std::copy(mat.data_, mat.data_ + N, data_);
}

Matrix::~Matrix() {
    delete[] data_;
}

double& Matrix::operator()(int i, int j) {
    if (i >= m_ || j >= n_)
        throw std::out_of_range("Indices out of bounds.");
    return data_[i*n_ + j];
}

double Matrix::operator()(int i, int j) const {
    if (i >= m_ || j >= n_)
        throw std::out_of_range("Indices out of bounds.");
    return data_[i*n_ + j];
}

Matrix& Matrix::zeros() {
    for (int i = 0; i < m_*n_; i++) {
        data_[i] = 0.;
    }
    return *this;
}
// Matrix DEFINITION END

// Cube DEFINITION START
Cube::Cube(int m,  int n,  int o)
        : m_ (m)
        , n_ (n)
        , o_ (o)
{
    if (m == 0 || n == 0 || o == 0)
        throw std::out_of_range("m, n, and o should be > 0");
    data_ = new double[m * n * o];
}

Cube::Cube(int m,  int n,  int o, double val) : Cube::Cube(m, n, o) {
    for (int i = 0; i < m*n*o; i++) {
        data_[i] = val;
    }
}

Cube::Cube(const Cube &cube) {
    m_ = cube.m_;
    n_ = cube.n_;
    o_ = cube.o_;
     int N = m_ * n_ * o_;
    data_ = new double[N];
    std::copy(cube.data_, cube.data_ + N, data_);
}

Cube::~Cube() {
    delete[] data_;
}

double& Cube::operator()(int i,  int j,  int k) {
    if (i >= m_ || j >= n_ || k >= o_)
        throw std::out_of_range("Indices out of bounds.");
    return data_[i*n_*o_ + j*o_ + k];
}

double Cube::operator()(int i,  int j,  int k) const {
    if (i >= m_ || j >= n_ || k >= o_)
        throw std::out_of_range("Indices out of bounds.");
    return data_[i*n_*o_ + j*o_ + k];
}

Cube& Cube::zeros() {
    for (int i = 0; i < m_*n_*o_; i++) {
        data_[i] = 0.;
    }
    return *this;
}
// Cube DEFINITION END
