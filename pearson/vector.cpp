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
    for (auto i{0}; i < size; i+=4)
    {
        data[i] = other.data[i];
        data[i+1] = other.data[i+1];
        data[i+2] = other.data[i+2];
        data[i+3] = other.data[i+3];
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
    double sum{0};

    for (auto i{0}; i < size; i+=4)
    {
        sum += data[i];
        sum += data[i+1];
        sum += data[i+2];
        sum += data[i+3];
    }

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
    __m256d sum = _mm256_setzero_pd();

    for(int i = 0; i < size; i+=4)
    {
        __m256d a = _mm256_loadu_pd(&data[i]);
        __m256d b = _mm256_loadu_pd(&rhs[i]);

        __m256d product = _mm256_mul_pd(a,b);

        sum = _mm256_add_pd(sum,product);
    }

    double result[4];
    _mm256_storeu_pd(result,sum);
    return result[0] + result[1] + result[2] + result[3];
}