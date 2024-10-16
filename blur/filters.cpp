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
        auto x_size {dst.get_x_size()};
        auto y_size {dst.get_y_size()};
        double w[Gauss::max_radius]{};
        Gauss::get_weights(radius, w); // depends on the radius not on the x or y

        for (auto x{0}; x < x_size; x++)
        {
            for (auto y{0}; y < y_size; y++)
            {
                auto r_x{w[0] * dst.r(x, y)}, g_x{w[0] * dst.g(x, y)}, b_x{w[0] * dst.b(x, y)}, n_x{w[0]};
                auto r_y{w[0] * scratch.r(x, y)}, g_y{w[0] * scratch.g(x, y)}, b_y{w[0] * scratch.b(x, y)}, n_y{w[0]};

                for (auto wi{1}; wi <= radius; wi++)
                {
                    auto wc{w[wi]};
                    auto x2{x - wi};
                    auto y2{y - wi};

                    if (x2 >= 0)
                    {
                        r_x += wc * dst.r(x2, y);
                        g_x += wc * dst.g(x2, y);
                        b_x += wc * dst.b(x2, y);
                        n_x += wc;
                    }
                    x2 = x + wi;
                    if (x2 < x_size)
                    {
                        r_x += wc * dst.r(x2, y);
                        g_x += wc * dst.g(x2, y);
                        b_x += wc * dst.b(x2, y);
                        n_x += wc;
                    }
                    if (y2 >= 0)
                    {
                        r_y += wc * scratch.r(x, y2);
                        g_y += wc * scratch.g(x, y2);
                        b_y += wc * scratch.b(x, y2);
                        n_y += wc;
                    }
                    y2 = y + wi;
                    if (y2 < dst.get_y_size())
                    {
                        r_y += wc * scratch.r(x, y2);
                        g_y += wc * scratch.g(x, y2);
                        b_y += wc * scratch.b(x, y2);
                        n_y += wc;
                    }
                }
                scratch.r(x, y) = r_x / n_x;
                scratch.g(x, y) = g_x / n_x;
                scratch.b(x, y) = b_x / n_x;

                dst.r(x, y) = r_y / n_y;
                dst.g(x, y) = g_y / n_y;
                dst.b(x, y) = b_y / n_y;
            }
        }

        return dst;
    }

}
