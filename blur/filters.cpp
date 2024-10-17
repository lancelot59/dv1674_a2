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

        // 2 Gauss operation before the loop not inside
        double w[Gauss::max_radius]{};
        Gauss::get_weights(radius, w);
        // 3 move this outside of the loop
        auto w0 {w[0]};

        //4 Moved those calling for those arrays outside of the loop. They are called once instead of multiple times
        auto R_dst = dst.get_R(), G_dst = dst.get_G(), B_dst = dst.get_B();
        auto R_scratch = scratch.get_R(), G_scratch = scratch.get_G(), B_scratch = scratch.get_B();

        // dst loop
        for (auto x{0}; x < x_size; x++)
        {
            for (auto y{0}; y < y_size; y++)
            {
                auto index {y * x_size + x};
                auto r{w0 * R_dst[index]}, g{w0 * G_dst[index]}, b{w0 * B_dst[index]}, n{w0};

                for (auto wi{1}; wi <= radius; wi++)
                {
                    auto wc{w[wi]};
                    auto x2{x - wi};

                    if (x2 >= 0)
                    {
                        auto left_index = y * x_size + x2;
                        r += wc * R_dst[left_index];
                        g += wc * G_dst[left_index];
                        b += wc * B_dst[left_index];
                        n += wc;
                    }
                    x2 = x + wi;
                    if (x2 < x_size)
                    {
                        auto right_index = y * x_size + x2;
                        r += wc * R_dst[right_index];
                        g += wc * G_dst[right_index];
                        b += wc * B_dst[right_index];
                        n += wc;
                    }
                }
                // those one stays because R, G, B are const
                scratch.r(x, y) = r / n;
                scratch.g(x, y) = g / n;
                scratch.b(x, y) = b / n;
            }
        }

        for (auto x{0}; x < x_size; x++)
        {
            for (auto y{0}; y < y_size; y++)
            {
                auto index {y * x_size + x};
                auto r{w0 * R_scratch[index]}, g{w0 * G_scratch[index]}, b{w0 * B_scratch[index]}, n{w0};

                for (auto wi{1}; wi <= radius; wi++)
                {
                    auto wc{w[wi]};
                    auto y2{y - wi};

                    if (y2 >= 0)
                    {
                        auto left_index = y2 * x_size + x;
                        r += wc * R_scratch[left_index];
                        g += wc * G_scratch[left_index];
                        b += wc * B_scratch[left_index];
                        n += wc;
                    }
                    y2 = y + wi;
                    if (y2 < y_size)
                    {
                        auto right_index = y2 * x_size + x;
                        r += wc * R_scratch[right_index];
                        g += wc * G_scratch[right_index];
                        b += wc * B_scratch[right_index];
                        n += wc;
                    }
                }
                dst.r(x, y) = r / n;
                dst.g(x, y) = g / n;
                dst.b(x, y) = b / n;
            }
        }

        return dst;
    }

}
