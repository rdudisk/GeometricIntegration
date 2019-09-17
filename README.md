# Geometric Integration

Geomi is a c++ library for numerical integration of physical models preserving the underlying geometrical structures.
This library is header only and heavily relies on templates.
Internal representation of mathematical objects make use of the
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) C++ library,
and non-linear systems solvers are based on the NOX package from the
[Trilinos](https://trilinos.github.io/) project.

## Variational methods

For Lagrangian systems, the module `Variational` implements the paradygm for variational symplectic methods
and a number of single step and multiple steps methods.
The case of a Lie group action invariant configuration manifold is not covered yet.
An [example](https://github.com/rdudisk/GeometricIntegration/tree/master/examples/Kepler)
is given for the two-bodies problem.

## Runge-Kutta Muthe-Kaas

The `RKMK` module contains classes implementing the Runge-Kutta Munthe-Kaas method principle.
Here is a concrete [example](https://github.com/rdudisk/GeometricIntegration/tree/master/examples/RigidBody)
for the rigid body problem.

## Install

Hints about installation on debian-based systems can be found
[here](https://github.com/rdudisk/GeometricIntegration/blob/master/install-process.txt).
