# es-flow

Atmospheric and marine flow characterisation tools.

A C++ port of Ocean Array Systems' original MATLAB-based flow analysis tools. This will allow rapid processing and fewer obstacles to deployment.


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


## First party dependencies

Other OAS-owned repsitories which are used by es-flow:

**es-instruments-x** Provides ascii and binary file readers for specific instrument types, from given manufacturers.

**utils** Provides utilities for general data processing and management.

**tg-engine** Provides numerical application of the biot savart equation for a collection of line vortices.


## Third party dependencies
 
We are currently using:
 
[**Intel MKL**]() to provide FFT and other performance primitives.
 
[**ceres-solver**](http://ceres-solver.org/index.html#), a well supported project by Google, is used for nonlinear least squares optimisation.

[**Eigen**](http://eigen.tuxfamily.org/) provides a linear algebra library. It isn't as consistent with MATLAB's API as armadillo, but is used extensively in ceres-solver so might be sensible to use in other areas of development too.

[**matio**](https://github.com/tbeu/matio) read and write tools for MATLAB .mat format files, including recent v7.3 (HDFS) file formats. Much higher level than writing the HDF5 files ourselves.

[**cxxopts**](https://github.com/jarro2783/cxxopts) argument parser for C++11 under the MIT license (NB most "standard" parsers are under GNU!!!).
 
[**glog**](https://github.com/google/glog) google's asynchronous logging library, used for logging to file.
 
[**bib2reSTcitation**](https://github.com/cykustcc/bib2reSTcitation) is a tool for converting .bib files to .txt files formatted as ReStructured Text references (useful for document generation compatible with sphinx).

[**sphinx**](http://www.sphinx-doc.org/en/1.5.1/) is the document generator; it'll take the .rst documentation files and turn them into formatted documentation, in either latex or HTML form.

[**sphinx_rtd_theme**](https://github.com/snide/sphinx_rtd_theme) gives us the excellent looking ReadTheDocs theme for our HTML documentation.

We're not yet committed to any of the following, but a range of possibly useful libraries is::
 
[**Armadillo**](http://arma.sourceforge.net) provides an extensive linear algebra library with a MATLAB-like API.

[**CppNumericalSolvers**](https://github.com/PatWie/CppNumericalSolvers) provides a directly analagous alternative to MATLAB's `fminsearch()`.

[**Linterp**](http://rncarpio.github.io/linterp/) provides a interpolation of gridded and unstructured data in N dimensions.

[**Tino Kluge**](http://kluge.in-chemnitz.de/opensource/spline/) maintains a spline interpolant library with linear extrapolation.

[**eigen-matio**](https://github.com/tesch1/eigen-matio) could be useful for reading and writing eigen matrix types to mat files.

### Third party library installation (OSX)

**Intel MKL:**
Download the Intel MKL library packages. Click on the icon and follow installation instructions. You'll need the administrator password. The tools are installed in `/opt/intel/`.
The `include` directory is `/opt/intel/include`.

**matio:**
Whatever you do, don't try to fork and build from source - the autoconf is complex and not suitable for OSX. Luckily there's a brew formula:
```bash
brew install homebrew/science/libmatio --with-hdf5
```
**ceres-solver including eigen and glog dependencies:**
```bash
brew install homebrew/science/ceres-solver
```
**cxxopts:**
Not necessary to install if simply deploying executables, as it's a header only library. To build es-flow, cxxopts must be installed alongside es-flow. From the es-flow root directory:
```bash
cd ..
git clone https://github.com/jarro2783/cxxopts
```
Then using cmake to build es-flow will find the headers correctly.

**Sphinx and sphinx_rtd_theme**
```bash
brew install python
pip install Sphinx
pip install sphinx_rtd_theme
```

### Third party library installation (Linux)


## Compilation

A cross-platform compilation file is provided using cmake.

Build process includes MATLAB based mex files for library functionality - requiring MATLAB to be installed on the build machine for linking purposes.


## Application Configuration

### Configuration

A 'default' configuration (app/config/default.py) is included which lists all settings and descriptions.

NOTE: Private access passwords and keys SHOULD NOT BE STORED IN A CODE REPOSITORY. Set your server environment up to define these keys (and other custom settings) via environment variables.

### Environment Variables

The following environment variables are required for the app to operate:

| **Name** | **Example** | **Description** |
| --- | --- | --- |
| TEST_DATA_DIR| /Users/thc29/Source/OAS/es-flow/test/test_data | The file path (absolute or relative to the application working directory) of the test data directory.|
### Requirements

The `requirements.txt` file is used to install the app's package dependencies. It is recommended that the versions and packages specified are updated in this file, as upgrades to packages should be considered part of version control, and be unit tested prior to deploymment.

The deployment process (below) describes the use of pip to ensure that the packages/versions actually installed match those in `requirements.txt`.

## Unit Testing

Unit testing is done with the google test framework and requires the following data files/folders to be installed in the TEST_DATA_DIR directory (see 'Environment Variables'):

- **fake_lidar_basic.mat** (generated by script `make_fake_lidar_basic.m`) which contains example lidar data for a basic lidar data structure (defined within the script). This .mat file is large (assumes high resolution temporal data over a long time). Actual flow profiles are not really valid; it's mostly just randomised data for unit testing purposes."

- **


