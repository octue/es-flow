=================
EnvironmentSTUDIO
=================

EnvironmentSTUDIO flow (**es-flow** for short) is a library for analysing and modeling atmospheric and marine boundary layers.

While **es-flow** can be used for any turbulent boundary layer analysis, its main focus is for
wind and tidal site characterisation - generating the 'best fit' of analytical models to measured velocity data.

Aims
====

The goal is to provide: [#f1]_

a. Profiles (for given input parameters) of
    - Mean Velocity (using power law, logarithmic law, MOST or Lewkowicz relations)
    - Mean Veer (Using Monin-Obukhov approach)
    - Reynolds Stress ``u'w'`` (using Lewkowicz relations)
    - Reynolds Stress ``u'u'``, ``u'v'``, ``u'w'``, ``v'v'``, ``v'w'``, ``w'w'`` (using the :doc:`adem`)
    - Spectra ``Sij`` (using Kaimal, von Karman)
    - Spectra ``Sij`` (using Attached-Detached Eddy Method)
    - Integral turbulent intensity and lengthscale ``I``, ``l``.

b. Best fit parameter sets (
    - given measured profiles from an instrument), to describe the above quantities analytically.

In future, generation of artificial flow fields for simulation purposes (eg ``.wnd`` fields for FVM/panel/BEM models,
inlets to DES, etc.) might be considered.

We use the `GitHub Issue Tracker <https://github.com/octue/es-flow>`_ to manage bug reports and feature requests.


Users
=====

At `Octue <http://www.octue.com>`_, **es-flow** is used to:

  * Provide a basis for developing and validating new processes, like the :doc:`adem`, for characterising Atmospheric Boundary Layers.
  * Process LiDAR and ADCP datasets from raw instrument files.
  * Apply windowed analyses to time series data, for flow characterisation.
  * Generate load cases for FAST, Bladed and our own TurbineGRID aerodynamic wind turbine analysis tools.

We'd like to hear about your use case. Please get in touch!


.. toctree::
   :maxdepth: 1
   :hidden:

   self
   adem
   examples
   file_formats
   installation
   license
   version_history
   bibliography
   library_api/library_root



.. rubric:: Footnotes

.. [#f1] Not all of this functionality is implemented yet!
