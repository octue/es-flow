.. _chapter-installation:

============
Installation
============


Third party library installation
================================

Intel MKL
---------

.. ATTENTION::
   If you don't wish to use Intel MKL, or need to build for a non-intel architecture, please contact Octue.

Download the Intel MKL library packages. Click on the icon and follow installation instructions. You'll need the administrator password.

The tools are installed in ``/opt/intel/``, the ``include`` directory is ``/opt/intel/include``.


matio
-----

.. tabs::

   .. group-tab:: Mac OSX

      Whatever you do, don't try to fork and build from source - the autoconf is complex and not suitable for OSX. Luckily there's a brew formula:

      .. code-block::

        brew install homebrew/science/libmatio --with-hdf5

   .. group-tab:: Linux

      Please contact Octue for Linux installation instructions.

   .. group-tab:: Windows

      Please contact Octue for Windows installation instructions.


ceres-solver, eigen and glog
----------------------------

.. tabs::

   .. group-tab:: Mac OSX

      Google's ceres-solver also depends on glog and eigen, so we get three for the price of one:

      .. code-block::

        brew install homebrew/science/ceres-solver


   .. group-tab:: Linux

      Please contact Octue for Linux installation instructions.

   .. group-tab:: Windows

      Please contact Octue for Windows installation instructions.


Third party build requirements
==============================

.. ATTENTION::
    These dependencies are only required if you're building **es-flow** from source.


cxxopts
-------

.. tabs::

   .. group-tab:: Mac OSX

      To build **es-flow**, ``cxxopts`` must be placed alongside **es-flow**. From the **es-flow** root directory:

      .. code-block::

        cd ../thirdparty
        git clone https://github.com/jarro2783/cxxopts

        Then using cmake to build **es-flow** will find the headers correctly.

   .. group-tab:: Linux

      Please contact Octue for Linux installation instructions.

   .. group-tab:: Windows

      Please contact Octue for Windows installation instructions.


NumericalIntegration
--------------------

.. tabs::

   .. group-tab:: Mac OSX

        To build **es-flow**, NumericalIntegration must be placed alongside **es-flow**. From the **es-flow** root directory:

      .. code-block::

        cd ../thirdparty
        git clone https://github.com/thclark/NumericalIntegration

        Then using cmake to build **es-flow** will find the headers correctly.

   .. group-tab:: Linux

      Please contact Octue for Linux installation instructions.

   .. group-tab:: Windows

      Please contact Octue for Windows installation instructions.
