# es-flow

Atmospheric and marine flow characterisation tools.

A C++ port of Ocean Array Systems' original MATLAB-based flow analysis tools. This will allow rapid processing and fewer obstacles to deployment.

___
## Project structure and Code Style

The folder structure is arranged as:
```
./         Makefile, readme, licensing and configure scripts.
./src      General sources
./include  Header files that expose the public interface and are to be installed
./lib      Library build directory
./bin      Tools and examples build directory
./docs	   Documentation, manually and autogenerated
./test     Test suites that should be run during a `make test`
./matlab   The deprecated MATLAB based toolset
```

Code style, includes and project structure should conform to the [Google C++ style guide](https://google.github.io/styleguide/cppguide.html) 

___
## First party dependencies

Other OAS-owned repsitories which are used by es-flow:

**es-instruments-x** Provides ascii and binary file readers for specific instrument types, from given manufacturers.

**utils** Provides utilities for general data processing and management.

**tg-engine** Provides numerical application of the biot savart equation for a collection of line vortices.

___
## Third party dependencies
 
[**Armadillo**](http://arma.sourceforge.net) provides an extensive linear algebra library with a MATLAB-like API.

[**Eigen**](http://eigen.tuxfamily.org/) provides a linear algebra library. It isn't as consistent with MATLAB's API as armadillo, but is used extensively in the ceres examples so could be easier for later stage development.

[**ceres-solver**](http://ceres-solver.org/index.html#), a well supported project by Google, is used for nonlinear least squares optimisation.

[**CppNumericalSolvers**](https://github.com/PatWie/CppNumericalSolvers) provides a directly analagous alternative to MATLAB's `fminsearch()`.

### Third party library installation (OSX)

### Third party library installation (Linux)

___
## Compilation

A cross-platform compilation file is provided using cmake.

Build process includes MATLAB based mex files for library functionality - requiring MATLAB to be installed on the build machine for linking purposes.
