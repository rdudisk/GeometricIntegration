Geomi: Geometric Integration										 {#mainpage}
============================

Presentation
------------

Geomi is a `C++` library to perform numerical integration of dynamical systems
while preserving the underlying geometrical structures.
Geomi implements a number of different integrators in several independent modules.
It is header only and heavily relies on templates.
Internal representation of mathematical objects make use of the
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) C++ library,
and non-linear systems solvers are based on the NOX package from the
[Trilinos](https://trilinos.github.io/) project.

The main modules are:
- `RKMK`: Runge-Kutta Munthe-Kaas methods
- `Variational`:
  implements variational and Lie-variational methods for Lagrangian systems.

The design concepts of Geomi are the following:
- Ease of use. ...
  However, the library is customizable; in particular, developping new
  integrators within the framework is relatively easy.
- Performance. The choice of `C++` and template classes is aimed towards
  efficiency and speed.
  Although some of the integrators could be faster if implemented for specific
  dynamical systems, we think that despite the tradeoff in flexibility over
  performance this implementation still compares pretty well (still has to be
  proven though...)
- Modularity
- Code factorisation

Installation
------------

Hints about installation on debian-based systems can be found
[here](https://github.com/rdudisk/GeometricIntegration/blob/master/install-process.txt).

Examples
--------

An [example](https://github.com/rdudisk/GeometricIntegration/tree/master/examples/Kepler)
is given for the two-bodies problem.
Here is a concrete [example](https://github.com/rdudisk/GeometricIntegration/tree/master/examples/RigidBody)
for the rigid body problem.

Quick tutorial
--------------
