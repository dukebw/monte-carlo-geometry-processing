# Monte Carlo Geometry Processing: A Grid-Free Approach to PDE-Based Methods on Volumetric Domains

I implemented a bare-bones version of the core idea from "Monte Carlo Geometry
Processing: A Grid-Free Approach to PDE-Based Methods on Volumetric Domains" by
Sawhney and Crane.
Suppose we want to solve a PDE like the Laplace or Poisson equation on a geometric structure.
Sawhney and Crane’s idea is to avoid meshing challenging geometric structures with finite-element methods (FEM).
FEM can be time-consuming for detailed objects such as the fly on the left (below).
Instead, they use Monte Carlo integration to find an approximate solution.
The Monte Carlo approximation is noisy at first, but can be refined by taking more samples, as shown on the right (below).

![](images/mcgp-1-Figure1-1.png)
(Figure 1 Sawhney & Crane 2020)


## Background

In my implementation I focused on solving the Laplace equation.

<p align="center"><img src="https://rawgit.com/dukebw/monte-carlo-geometry-processing/main/svgs/bd19e4aa14e8e7729bd6f7be2cb9ec6e.svg?invert_in_darkmode" align=middle width=172.05074864999997pt height=14.611878599999999pt/></p>

To solve the Laplace equation, Sawhney and Crane use a technique called "walk on spheres".
The idea uses the fact that the solution of a harmonic function at point p is
the expected value of the first point on the boundary reached by a random walk
starting at p.
By symmetry, we can find this boundary point by sampling uniformly on spheres.
We start at x0, then successively sample boundary point until we get to `x_k`,
which is within some threshold <img src="./svgs/e545a3e80cb1165111423624206de470.svg?invert_in_darkmode" align=middle width=26.99150024999999pt height=22.831056599999986pt/> of the boundary (illustrated below).

![](images/mcgp-3-Figure3-1.png)
(Figure 3 Sawhney & Crane 2020)


## Code

I implemented two versions of walk on spheres: a 2D version `walk_on_circles` and a 3D version `walk_on_spheres`.
I put their respective implementations in `walk_on_circles.h/cpp` and `walk_on_spheres.h/cpp`.
Install the code with cmake as usual.
Then, run the `walk_on_circles` code from `main.cpp` with the command line

```
./monte-carlo-geometry-processing --nvertices 256 --max-steps 16 --num-walks 8
```.

`nvertices` controls the number of vertices in the planar mesh, `max-steps` is
the maximum number of walk on circles steps, and `num-walks` is the number of
walks on circles per vertex.
