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


Intel TBB
---------

.. ATTENTION::
   If you don't wish to use Intel TBB, or need to build for a non-intel architecture, please contact Octue.


.. tabs::

   .. group-tab:: Mac OSX

      TBB is installable via brew.

      .. code-block:: bash

          brew install tbb

   .. group-tab:: Linux

        Download the Intel TBB library packages. Click on the icon and follow installation instructions, ensuring that
        the TBBROOT environment variable is set. You'll need the administrator password.

        The tools are installed in ``/opt/intel/``, the ``include`` directory is ``/opt/intel/include``.

   .. group-tab:: Windows

        Follow the Intel TBB instructions, ensuring that the ``TBBROOT`` environment variable is set.

matio
-----

.. tabs::

   .. group-tab:: Mac OSX

      Building from source is possible using instructions on the matio home page, although unwieldy. The
      `next major release of matio <https://github.com/tbeu/matio/issues/133>`_ should bring CMake to the party,
      whereupon we'll add that to the build
      system transparently (downloading and building it if not found).

      In the meantime, you'll most likely have success with a brew formula, although `check this issue<>`_ if you
      experience runtime problems with HDF5 loading the signature files.

      .. code-block:: bash

          brew install libmatio

   .. group-tab:: Linux

      Please contact Octue for Linux installation help.

   .. group-tab:: Windows

      Please contact Octue for Windows installation help.


ceres-solver, eigen and glog
----------------------------

.. tabs::

   .. group-tab:: Mac OSX

      Google's ceres-solver also depends on glog and eigen, so we get three for the price of one.

      .. code-block:: bash

          brew install homebrew/science/ceres-solver


   .. group-tab:: Linux

      Please contact Octue for Linux installation help.

   .. group-tab:: Windows

      Please contact Octue for Windows installation help.


Third party build requirements
==============================

.. ATTENTION::
    These dependencies are only required if you're building **es-flow** from source.


cxxopts
-------

.. tabs::

   .. group-tab:: Mac OSX

      To build **es-flow**, ``cxxopts`` must be placed alongside **es-flow**. From the **es-flow** root directory:

      .. code-block:: bash

          cd ../thirdparty
          git clone https://github.com/jarro2783/cxxopts

      Then using cmake to build **es-flow** will find the headers correctly.

   .. group-tab:: Linux

      Please contact Octue for Linux installation help.

   .. group-tab:: Windows

      Please contact Octue for Windows installation help.


NumericalIntegration
--------------------

.. tabs::

   .. group-tab:: Mac OSX

      To build **es-flow**, NumericalIntegration must be placed alongside **es-flow**. From the **es-flow** root directory:

      .. code-block:: bash

          cd ../thirdparty
          git clone https://github.com/thclark/NumericalIntegration

      Then using cmake to build **es-flow** will find the headers correctly.

   .. group-tab:: Linux

      Please contact Octue for Linux installation help.

   .. group-tab:: Windows

      Please contact Octue for Windows installation help.
