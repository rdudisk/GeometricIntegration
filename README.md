Geomi: Geometric Integration										 {#mainpage}
============================

## Presentation

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

## Installation

In order for Geomi to compile and run, some dependencies have to be installed.
Hints about the process are given for debian-based systems.

### NOX

NOX is a package from the [Trilinos](https://trilinos.github.io/) project.
To install Trilinos, run

    # apt-get install trilinos-all-dev trilinos-doc

### Eigen 3

Eigen is a C++ library for linear algebra.
Information about its installation can be found on the
[official website](http://eigen.tuxfamily.org/index.php?title=Main_Page).
For debian users, there is a package available

    # apt-get install libeigen3-dev libeigen3-doc

or, for the up-to-date version, the official GitHub

    $ cd ~
    $ git clone https://github.com/eigenteam/eigen-git-mirror
    $ cmake eigen-git-mirror
    # make install

If includes directives do not work when you compile, you probably need to
create a link (you need to adapt depending on the actual installation directory)

    # ln -s /usr/include/eigen3/Eigen /usr/local/include/Eigen

<!--
	Nécessaire ? Essayer sans d'abord
	# apt-get install libmrpt-dev
-->

### Doxygen

If you want to be able to compile the documentation, you need to install
[Doxygen](http://doxygen.nl/)

    # apt-get install doxygen*

## Examples

For the `Variational` module, an
[example](https://github.com/rdudisk/GeometricIntegration/tree/master/examples/Kepler)
is given for the two-bodies problem.
Regarding the `RKMK` module, a solution for the rigid body problem is given
[here](https://github.com/rdudisk/GeometricIntegration/tree/master/examples/RigidBody/RKMK).

## Quick tutorial

An introduction on how to use the library can be found on this
[page](https://rdudisk.github.io/GeometricIntegration/doc/html/quick_tutorial.html).

## Documentation

Online documentation can be found
[here](https://rdudisk.github.io/GeometricIntegration/doc/html/index.html).
It can also be build from the source by running `$ make doxygen`

