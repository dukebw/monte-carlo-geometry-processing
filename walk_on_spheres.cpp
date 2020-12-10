#include "walk_on_spheres.h"
#include "Eigen/Core"
#include "igl/parallel_for.h"
#include "igl/point_mesh_squared_distance.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <random>

static float
walk_on_spheres_from_point(const Eigen::MatrixXf& V,
                           const Eigen::MatrixXi& boundary_edges,
                           int num_walks,
                           int max_steps,
                           int Vidx)
{

        float stopping_tolerance = 0.05f;

        std::random_device rd;
        std::mt19937 gen{ rd() };
        std::normal_distribution<float> dis{ 0.0f, 1.0f };

        float accum = 0.0f;
        int num_successful_walks = 0;
        Eigen::Vector3f x0 = V.row(Vidx);
        for (int i = 0; i < num_walks; i++) {
                Eigen::MatrixXf x{ 1, 3 };
                x.row(0) = x0;

                for (int steps = 0; steps < max_steps; ++steps) {
                        float radius = std::numeric_limits<float>::max();

                        Eigen::VectorXf squared_distances;
                        Eigen::VectorXi closest_edge;
                        Eigen::MatrixXf closest_points;
                        igl::point_mesh_squared_distance(x,
                                                         V,
                                                         boundary_edges,
                                                         squared_distances,
                                                         closest_edge,
                                                         closest_points);

                        radius = sqrtf(squared_distances(0));
                        if (radius < stopping_tolerance) {
                                accum += closest_points(0, 2);
                                ++num_successful_walks;
                                break;
                        }

                        /**
                         * NOTE(brendan): sample from a three-dimensional Gaussian.
                         * This is equivalent to sampling a 3-vector with a random
                         * orientation in space.
                         */
                        Eigen::Vector3f dx{ dis(gen), dis(gen), dis(gen) };
                        dx /= dx.norm() + 0.0001f;
                        x.row(0) += radius * dx;
                }
        }

        if (num_successful_walks == 0)
                return 0.0f;

        return accum / num_successful_walks;
}

void
walk_on_spheres(Eigen::VectorXf& Z,
                const Eigen::MatrixXf& V,
                const Eigen::MatrixXi& boundary_edges,
                int num_walks,
                int max_steps)
{
        const auto walk_on_spheres_parallel = [&](const int Vidx) {
                Z(Vidx) = walk_on_spheres_from_point(
                  V, boundary_edges, num_walks, max_steps, Vidx);
        };
        igl::parallel_for(V.rows(), walk_on_spheres_parallel, 1000);
}
