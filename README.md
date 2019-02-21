# es-flow

Atmospheric and marine flow characterisation tools, implementing the 'Attached-Detached Eddy Method' for boundary layer analysis.


## Project structure and Code Style

The folder structure is arranged as:

- `./`         CMake file, readme, licensing, requirements (for python) and git configuration files
- `./docs`	   Documentation in `.rst` form (including autogenerated library docs)
- `./matlab`   The deprecated MATLAB based toolset
- `./scripts`  Supporting scripts (eg `make_docs.py` to autogenerate the library API docs using sphinx and breathe)
- `./source`   Source code
- `./test`     Unit tests

CMake builds in its own directories. However, once built, files may be moved / copied to the following locations for inclusion in other projects:
- `./include`  Header files that expose the public interface to the library
- `./lib`      Library files (`.a`, `.so`, `.dylib`, `.dll` etc)
- `./bin`      Main executable, plus tools and examples


Code style, includes and project structure should conform to the [Google C++ style guide](https://google.github.io/styleguide/cppguide.html) 


## Release Flow

Steps to create a new release:

- Create a new branch, from master, entitled `release<version-tag>` and push to remote. Version tags are semantic:
```
git checkout master
git checkout -b release-v0.1.1-alpha.1
git push --set-upstream origin release-v0.1.1-alpha.1
```
- TravisCI automatically creates a tagged draft release on github. In the above example, the tag is `v0.1.1-alpha.1`.
- Any commits into, or merges to, this branch, will trigger automatic builds on TravisCI. The documentation gets built, binaries compiled, tests run. Assets (e.g. precompiled binaries and HTML documentation) are assembled into a `zip` file and appended to the release. These are uploaded for each commit. 
- Any new features, fixes etc on your development branch should first be merged into the release, not into master. 
- Once the release is prepared, checked, and all builds have completed, publish the release from GitHub.
- Once the release is published, the feature branch may finally be merged into master.


## Third party dependencies

### Currently in use

See documentation
 
[**Intel MKL**]() to provide FFT and other performance primitives, enhancing Eigen.
 
[**ceres-solver**](http://ceres-solver.org/index.html#), a well supported project by Google, is used for nonlinear least squares optimisation.

[**Eigen**](http://eigen.tuxfamily.org/) provides a linear algebra library. It isn't as consistent with MATLAB's API as armadillo, but is used extensively in ceres-solver and elsewhere so is our selected library of choice.

[**NumericalIntegration**](https://github.com/thclark/NumericalIntegration) NB This fork avoids compiling against the GPL licensed MPFRC++ library. Provides a numerical integration module for eigen. It's weird that it's not an unsupported module in eigen ([see this thread](https://bitbucket.org/eigen/eigen/pull-requests/109/numerical-integration-module-for-eigen/diff)) given that NumericalDiff *is* part of eigen.

[**matio**](https://github.com/tbeu/matio) read and write tools for MATLAB .mat format files, including recent v7.3 (HDFS) file formats. Much higher level than writing the HDF5 files ourselves.

[**cxxopts**](https://github.com/jarro2783/cxxopts) argument parser for C++11 under the MIT license (NB most "standard" parsers are under GNU!!!).
 
[**glog**](https://github.com/google/glog) google's asynchronous logging library, used for logging to file.
 
[**sphinx**](http://www.sphinx-doc.org/en/1.5.1/) is the document generator; it'll take the .rst documentation files and turn them into formatted documentation, in either latex or HTML form.

[**sphinx_rtd_theme**](https://github.com/snide/sphinx_rtd_theme) gives us the excellent looking ReadTheDocs theme for our HTML documentation.



### For consideration and possible future use

We're not yet committed to any of the following, but a range of possibly useful libraries is::

[**gflags**](https://github.com/gflags/gflags) is an optional alternative to cxxopts, which is more widely supported and already used by ceres (i.e. installed already as part of the build chain). We have an issue open to change over (see #18) to gflags.

[**CppNumericalSolvers**](https://github.com/PatWie/CppNumericalSolvers) provides a directly analagous alternative to MATLAB's `fminsearch()`.

[**Linterp**](http://rncarpio.github.io/linterp/) provides a interpolation of gridded and unstructured data in N dimensions.

[**Tino Kluge**](http://kluge.in-chemnitz.de/opensource/spline/) maintains a spline interpolant library with linear extrapolation.

[**eigen-matio**](https://github.com/tesch1/eigen-matio) could be useful for reading and writing eigen matrix types to mat files.


## Compilation

A cross-platform compilation file is provided using cmake.

### MATLAB MEX
**DEPRECATION WARNING:** The MathWorks really make it extremely difficult to support integrations. *MATLAB mex files are almost completely platform- and matlab-version- dependent. So you basically have to recompile for every MATLAB release, on every platform, clearly a nightmare. MATLAB does however allow you to invoke C library functions directly (see `loadlibrary` in the MATLAB help) so we may consider writing a header exposing C-style library API for use with MATLAB, to allow us to continue supporting MATLAB in a practical way.*

The CMake build process includes MATLAB based mex files for library functionality - requiring MATLAB to be installed on the build machine for linking purposes. This is unsupported as of January 2019.


## Documentation

Documentation resides in the `./docs` directory, as `*.rst` files. The instruction manual includes an automatically generated subfolder of .rst files detailing the API of the library.

### Creating the bibliography

For sphinx to create a bibliography, the `bibliography.rst` file needs to contain, in RestructuredText format, the references used.

However, you may have  BibTeX, rather than .rst, references. To convert, you can use [**bib2reSTcitation**](https://github.com/cykustcc/bib2reSTcitation), a handy tool for converting `.bib` files to `.rst` files.

### Building documentation

The build steps are as follows:

- Do any conversions of reference formats that you need to, ensuring `bibliography.rst` is up to date with all references.
- Use `breathe` to parse the library and generate `doxygen` files (in `xml` form)
- Use `exhale` to convert doxygen `xml` files to `rst` documentation
- Use `sphinx` (with `mathjax` to convert AMS LaTeX to rendered equations and `sphinx-rtd-theme` for prettiness) to convert all the `rst` files to `html`.

### Build Environment

The build scripts are in python 2.7, because sphinx hasn't been fully moved to python 3 yet. Sigh!

If developing docs, see the `.travis.yml` file for the doc build steps - you should be able to repeat these or similar on your machine.

If you're using pyenv locally and have other versions of sphinx installed, you also need to add a shim to pyenv so that sphinx calls the right python installation. Copy `scripts/sphinx-build` to the pyenv shim directory, in my case its `/Users/thc29/.pyenv/shims/`, and edit the last two lines of it to your settings.


## Unit Testing

Unit testing is done with the google test framework and requires the following data files/folders to be installed in the `TEST_DATA_DIR` directory (see 'Environment Variables'):

- **fake_lidar_basic.mat** (generated by script `make_fake_lidar_basic.m`) which contains example lidar data for a basic lidar data structure (defined within the script). This .mat file is large (assumes high resolution temporal data over a long time). Actual flow profiles are not really valid; it's mostly just randomised data for unit testing purposes.


### Test Environment

The following environment variables are required for the app to operate:

| **Name** | **Example** | **Description** |
| --- | --- | --- |
| `TEST_DATA_DIR` | `/Users/thc29/Source/octue/es-flow/test/test_data` | The file path (absolute or relative to the application working directory) of the test data directory, WITHOUT a trailing slash.|
