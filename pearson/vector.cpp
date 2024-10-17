/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "vector.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <immintrin.h>

Vector::Vector()
    : size{0}, data{nullptr}
{
}

Vector::~Vector()
{
    if (data)
    {
        delete[] data;
    }

    size = 0;
}

Vector::Vector(unsigned size)
    : size{size}, data{new double[size]}
{
}

Vector::Vector(unsigned size, double *data)
    : size{size}, data{data}
{
}

Vector::Vector(const Vector &other)
    : Vector{other.size}
{
    for (int i = 0; i < size; i+=4)
    {
        __m256d dataVec = _mm256_loadu_pd(&other.data[i]);
        _mm256_storeu_pd(&data[i], dataVec);
    }
}

unsigned Vector::get_size() const
{
    return size;
}

double *Vector::get_data()
{
    return data;
}

double Vector::operator[](unsigned i) const
{
    return data[i];
}

double &Vector::operator[](unsigned i)
{
    return data[i];
}

double Vector::mean() const
{
    __m256d sumVec = _mm256_setzero_pd();
    int i = 0;

    for (; i <= size -4; i+=4)
    {
        __m256d dataVec = _mm256_loadu_pd(&data[i]);
        sumVec = _mm256_add_pd(sumVec, dataVec);
    }
    double sum = 0;
    double temp[4];
    _mm256_storeu_pd(temp, sumVec);
    sum += temp[0] + temp[1] + temp[2] + temp[3];
    return sum / static_cast<double>(size);
}

double Vector::magnitude() 
{
    auto dot_prod{dot(*this)};
    return std::sqrt(dot_prod);
}

Vector& Vector::operator/(double div)
{

    for (auto i{0}; i < size; i+=4)
    {
        data[i] /= div;
        data[i+1] /= div;
        data[i+2] /= div;
        data[i+3] /= div;
    }

    return *this;
}

Vector& Vector::operator-(double sub)
{

    for (auto i{0}; i < size; i+=4)
    {
        data[i] -= sub;
        data[i+1] -= sub;
        data[i+2] -= sub;
        data[i+3] -= sub;
    }

    return *this;
}

double Vector::dot(Vector& rhs) const
{
    double result{0};

    for (auto i{0}; i < size; i+=4)
    {
        result += data[i] * rhs[i];
        result += data[i+1] * rhs[i+1];
        result += data[i+2] * rhs[i+2];
        result += data[i+3] * rhs[i+3];
    }

    return result;
}