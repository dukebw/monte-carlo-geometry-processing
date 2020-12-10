#include "walk_on_circles.h"

#include "Eigen/Core"
#include "igl/PI.h"
#include "igl/parallel_for.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <random>

static float
closest_point_on_segment(const Eigen::MatrixXf& bs,
                         const Eigen::Vector2f& x,
                         int seg_idx,
                         float R)
{
        Eigen::Vector2f seg0{ bs(seg_idx, 0), bs(seg_idx, 1) };
        Eigen::Vector2f seg1{ bs(seg_idx, 2), bs(seg_idx, 3) };

        Eigen::Vector2f u = seg1 - seg0;
        float t = ((x - seg0).transpose() * u);
        t /= u.squaredNorm();
        t = std::clamp(t, 0.0f, 1.0f);

        Eigen::Vector2f p = ((1 - t) * seg0) + (t * seg1);
        float distance = (x - p).norm();

        return std::min(R, distance);
}

static float
walk_on_circles_for_point(const Eigen::Vector2f& x0,
                          const Eigen::MatrixXf& bndary_segs,
                          const std::function<float(const Eigen::Vector2f&)> g,
                          int num_walks,
                          int max_steps)
{
        float stopping_tolerance = 0.01f;

        std::random_device rd;
        std::mt19937 gen{ rd() };
        std::uniform_real_distribution<float> dis{ 0.0f, 2.0f * igl::PI };

        float accum = 0.0f;
        Eigen::Vector2f x;
        int num_successful_walks = 0;
        for (int i = 0; i < num_walks; i++) {
                Eigen::Vector2f x = x0;
                for (int steps = 0; steps < max_steps; ++steps) {
                        float radius = std::numeric_limits<float>::max();
                        for (int seg_idx = 0; seg_idx < bndary_segs.rows(); ++seg_idx) {
                                radius = closest_point_on_segment(
                                  bndary_segs, x, seg_idx, radius);
                        }

                        if (radius < stopping_tolerance) {
                                accum += g(x);
                                ++num_successful_walks;
                                break;
                        }

                        /**
                         * NOTE(brendan): sample a uniform random point on the circle
                         * centered at x.
                         * The radius of the circle is the distance of the closest point
                         * on the boundary to x.
                         */
                        float theta = dis(gen);
                        x(0) += radius * cosf(theta);
                        x(1) += radius * sinf(theta);
                }
        }

        if (num_successful_walks == 0)
                return 0.0f;

        return accum / num_successful_walks;
}

void
walk_on_circles(Eigen::VectorXf& Z,
                const Eigen::MatrixXf& V,
                const Eigen::MatrixXf& boundary_segments,
                const std::function<float(const Eigen::Vector2f&)> g,
                int num_walks,
                int max_steps)
{
        const auto walk_on_circles_parallel = [&](const int Vidx) {
                Z(Vidx) =
                  walk_on_circles_for_point(Eigen::Vector2f{ V(Vidx, 0), V(Vidx, 1) },
                                            boundary_segments,
                                            g,
                                            num_walks,
                                            max_steps);
        };
        igl::parallel_for(V.rows(), walk_on_circles_parallel, 1000);
}
