#ifndef WALK_ON_CIRCLES_H
#define WALK_ON_CIRCLES_H
#include <Eigen/Core>
#include <functional> // for std::function

/**
 * Solve Laplace equation using walk on spheres in 2D.
 * Give the boundary condition in piecewise segments (boundary_segments).
 * Pass a function g to evaluate on the boundary.
 *
 * Inputs:
 *      Z #V vector of per-vertex solution values
 *      V #V by 3 list of mesh vertex positions
 *      boundary_segments #boundary_segments by 4 list of (x0, y0, x1, y1) segments
 *      g boundary condition function
 *      num_walks number of random walks on spheres per vertex
 *      max_steps maximum steps per walk on spheres
 */
void
walk_on_circles(Eigen::VectorXf& Z,
                const Eigen::MatrixXf& V,
                const Eigen::MatrixXf& boundary_segments,
                const std::function<float(const Eigen::Vector2f&)> g,
                int num_walks,
                int max_steps);

#endif /* WALK_ON_CIRCLES_H */
