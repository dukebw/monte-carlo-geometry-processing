#define EIGEN_NO_DEBUG
#include "./walk_on_circles.h"
#include "./walk_on_spheres.h"
#include "Eigen/Sparse"
#include "cxxopts.hpp"
#include "igl/boundary_facets.h"
#include "igl/colon.h"
#include "igl/cotmatrix.h"
#include "igl/jet.h"
#include "igl/min_quad_with_fixed.h"
#include "igl/opengl/glfw/Viewer.h"
#include "igl/point_mesh_squared_distance.h"
#include "igl/readOFF.h"
#include "igl/setdiff.h"
#include "igl/slice.h"
#include "igl/slice_into.h"
#include "igl/unique.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>

static Eigen::MatrixXf V;
static Eigen::MatrixXi F;
static Eigen::VectorXf Z;
static Eigen::MatrixXf boundary_segments;
static Eigen::MatrixXi boundary_edges;
static uint32_t n_verts;
static int max_steps;
static int num_walks;
static uint32_t num_updates;
static bool use_camel;

/**
 * NOTE(brendan): Solve Laplace equation using walk on spheres in 3D.
 * Use the z coordinate of the boundary points as boundary condition.
 */
static float
walk_on_spheres(int Vidx)
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

static float
boundary_condition(const Eigen::Vector2f& x)
{
        return fmod((6.0f * x.array()).floor().sum(), 2.0f);
}

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

/**
 * NOTE(brendan): Solve Laplace equation using walk on spheres in 2D.
 * Give the boundary condition in piecewise segments (bndary_segs).
 * Pass a function g to evaluate on the boundary.
 */
static float
walk_on_circles(const Eigen::Vector2f& x0,
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

                        float theta = dis(gen);
                        x(0) += radius * cosf(theta);
                        x(1) += radius * sinf(theta);
                }
        }

        if (num_successful_walks == 0)
                return 0.0f;

        return accum / num_successful_walks;
}

static void
update_Z(void)
{
        if (use_camel) {
                const auto walk_on_spheres_parallel = [&](const int Vidx) {
                        Z(Vidx) = walk_on_spheres(Vidx);
                };
                igl::parallel_for(V.rows(), walk_on_spheres_parallel, 1000);
        } else {
                const auto walk_on_circles_parallel = [&](const int Vidx) {
                        Z(Vidx) =
                          walk_on_circles(Eigen::Vector2f{ V(Vidx, 0), V(Vidx, 1) },
                                          boundary_segments,
                                          boundary_condition,
                                          num_walks,
                                          max_steps);
                };
                igl::parallel_for(V.rows(), walk_on_circles_parallel, 1000);
        }
        ++num_updates;
}

static bool
key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
        if (key == '1') {
                Eigen::VectorXf Z_cma = Z;
                update_Z();
                Z = (((num_updates - 1) * Z_cma) + Z) / num_updates;
                viewer.data().set_data(
                  Z.cast<double>(), igl::COLOR_MAP_TYPE_TURBO, 256);
        } else if (key == '2') {
                num_updates = 0;
                update_Z();
                viewer.data().set_data(
                  Z.cast<double>(), igl::COLOR_MAP_TYPE_TURBO, 256);
        }
        return false;
}

int
main(int argc, char** argv)
{
        cxxopts::Options options("Monte Carlo Geometry Processing",
                                 "Monte Carlo Geometry Processing on a plane");

        options.add_options()                                                    //
          ("n,nvertices", "Num vertices (per edge)", cxxopts::value<uint32_t>()) //
          ("e,exit-before-gui",
           "Exit before GUI",
           cxxopts::value<bool>()->default_value("false")) //
          ("c,use-camel",
           "Use camel?",
           cxxopts::value<bool>()->default_value("false")) //
          ("m,max-steps",
           "Max walk on spheres steps",
           cxxopts::value<int>()->default_value("16")) //
          ("w,num-walks",
           "Number of walks on spheres",
           cxxopts::value<int>()->default_value("128")) //
          ;

        auto result = options.parse(argc, argv);

        n_verts = result["nvertices"].as<uint32_t>();
        bool should_exit_before_gui = result["exit-before-gui"].as<bool>();
        use_camel = result["use-camel"].as<bool>();
        max_steps = result["max-steps"].as<int>();
        num_walks = result["num-walks"].as<int>();

        uint32_t n_edges = n_verts - 1;
        float grid_spacing = 1.0 / n_verts;

        V.resize(n_verts * n_verts, 3);
        F.resize(2 * n_edges * n_edges, 3);

        if (use_camel) {
                igl::readOFF("../camelhead.off", V, F);
                igl::boundary_facets(F, boundary_edges);
        } else {
                /* NOTE(brendan): plane */
                for (int i = 0; i < n_verts; ++i) {
                        for (int j = 0; j < n_verts; ++j) {
                                uint32_t Vidx = (i * n_verts) + j;
                                V(Vidx, 0) = i * grid_spacing;
                                V(Vidx, 1) = j * grid_spacing;
                                V(Vidx, 2) = 0.0;
                        }
                }

                for (int i = 0; i < n_edges; ++i) {
                        for (int j = 0; j < n_edges; ++j) {
                                uint32_t Fidx = (i * 2 * n_edges) + (2 * j);
                                F(Fidx, 0) = (i * n_verts) + j;
                                F(Fidx, 1) = (i * n_verts) + (j + 1);
                                F(Fidx, 2) = ((i + 1) * n_verts) + (j + 1);

                                Fidx = (i * 2 * n_edges) + (2 * j) + 1;
                                F(Fidx, 0) = (i * n_verts) + j;
                                F(Fidx, 1) = ((i + 1) * n_verts) + (j + 1);
                                F(Fidx, 2) = ((i + 1) * n_verts) + j;
                        }
                }

                boundary_segments.resize(5, 4);
                boundary_segments << 0.5, 0.1, 0.9, 0.5, // x0, y0, x1, y1
                  0.5, 0.9, 0.1, 0.5,                    //
                  0.1, 0.5, 0.5, 0.1,                    //
                  0.5, 0.333, 0.5, 0.667,                //
                  0.333, 0.5, 0.667, 0.5;
        }

        Z.resize(V.rows());
        /* NOTE(brendan): alternative with monte carlo */
        update_Z();

        if (should_exit_before_gui)
                std::exit(EXIT_SUCCESS);

        // Plot the mesh with pseudocolors
        igl::opengl::glfw::Viewer viewer;

        viewer.callback_key_down = &key_down;

        viewer.data().set_mesh(V.cast<double>(), F);
        viewer.data().show_lines = false;
        viewer.data().set_data(Z.cast<double>(), igl::COLOR_MAP_TYPE_TURBO, 256);
        viewer.launch();
}
