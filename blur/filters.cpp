/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>

namespace Filter
{

    namespace Gauss
    {
        void get_weights(int n, double *weights_out)
        {
            for (auto i{0}; i <= n; i++)
            {
                double x{static_cast<double>(i) * max_x / n};
                weights_out[i] = exp(-x * x * pi);
            }
        }
    }

    Matrix blur(Matrix m, const int radius)
    {
        Matrix scratch{PPM::max_dimension};
        auto dst{m};
        // 1 adding getting date before the loop
        auto x_size {dst.get_x_size()};
        auto y_size {dst.get_y_size()};
        // 2 Gaus operation before the loop not inside
        double w[Gauss::max_radius]{};
        Gauss::get_weights(radius, w);
        // 3 move this outside of the loop
        auto w0 {w[0]};

        // 5 directly accessing R, G, and B - had to change the get method in MAtrix (deleted const)
        auto dst_R {dst.get_R()}, dst_B {dst.get_B()}, dst_G {dst.get_G()};
        auto scr_R {scratch.get_R()}, scr_B {scratch.get_B()}, scr_G {scratch.get_G()};

        // dst loop
        for (auto x{0}; x < x_size; x++)
        {
            for (auto y{0}; y < y_size; y++)
            {
                // 4 index is calculated before and then values are accessed by it (without calculations)
                auto index {x + y * x_size};
                auto r{w0 * dst_R[index]}, g{w0 * dst_G[index]}, b{w0 * dst_B[index]}, n{w0};

                for (auto wi{1}; wi <= radius; wi++)
                {
                    auto wc{w[wi]};
                    auto x2{x - wi};

                    if (x2 >= 0)
                    {
                        auto left_index {x2 + y * x_size};
                        r += wc * dst_R[left_index];
                        g += wc * dst_G[left_index];
                        b += wc * dst_B[left_index];
                        n += wc;
                    }
                    x2 = x + wi;
                    if (x2 < x_size)
                    {
                        auto right_index {x2 + y * x_size};
                        r += wc * dst_R[right_index];
                        g += wc * dst_G[right_index];
                        b += wc * dst_B[right_index];
                        n += wc;
                    }
                }
                scratch.r(index) = r / n;
                scratch.g(index) = g / n;
                scratch.b(index) = b / n;
            }
        }

        for (auto x{0}; x < x_size; x++)
        {
            for (auto y{0}; y < y_size; y++)
            {
                auto index = x + y * x_size;
                auto r{w0 * scr_R[index]}, g{w0 * scr_G[index]}, b{w0 * scr_B[index]}, n{w0};

                for (auto wi{1}; wi <= radius; wi++)
                {
                    auto wc{w[wi]};
                    auto y2{y - wi};

                    if (y2 >= 0)
                    {
                        auto left_index {x + y2 * x_size};
                        r += wc * scr_R[left_index];
                        g += wc * scr_G[left_index];
                        b += wc * scr_B[left_index];
                        n += wc;
                    }
                    y2 = y + wi;
                    if (y2 < y_size)
                    {
                        auto right_index {x + y2 * x_size};
                        r += wc * scr_R[right_index];
                        g += wc * scr_G[right_index];
                        b += wc * scr_B[right_index];
                        n += wc;
                    }
                }
                dst.r(index) = r / n;
                dst.g(index) = g / n;
                dst.b(index) = b / n;
            }
        }

        return dst;
    }

}