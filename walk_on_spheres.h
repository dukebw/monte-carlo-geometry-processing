#ifndef WALK_ON_SPHERES_H
#define WALK_ON_SPHERES_H
#include <Eigen/Core>

/**
 * Solve Laplace equation using walk on spheres in 3D.
 * Use the z coordinate of the boundary points as boundary condition.
 *
 * Inputs:
 *      Z #V vector of per-vertex solution values
 *      V #V by 3 list of mesh vertex positions
 *      boundary_edges #boundary_edges by 2 list of boundary edges for the mesh
 *      num_walks number of random walks on spheres per vertex
 *      max_steps maximum steps per walk on spheres
 */
void
walk_on_spheres(Eigen::VectorXf& Z,
                const Eigen::MatrixXf& V,
                const Eigen::MatrixXi& boundary_edges,
                int num_walks,
                int max_steps);

#endif /* WALK_ON_SPHERES_H */
