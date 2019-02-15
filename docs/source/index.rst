=================
EnvironmentSTUDIO
=================

EnvironmentSTUDIO flow (**es-flow** for short) is a library for analysing and modeling atmospheric and marine boundary layers.

While **es-flow** can be used for any turbulent boundary layer analysis, its main focus is for
wind and tidal site characterisation - generating the 'best fit' of analytical models to measured velocity data.

A key strength of **es-flow** is the :doc:`adem`. This extremely robust method allows users to:

    - determine detailed turbulence information from instruments like LiDAR
    - characterise turbulence and shear beyond tip height of even the biggest offshore wind turbines
    - characterise **coherent structure** in turbulence, crucially important for fatigue loading in wind turbines.


Aims
====

The **es-flow** library provides: [#f1]_

#. **Parameterised profiles**
    - Mean Velocity (using power law, logarithmic law, MOST or Lewkowicz relations)
    - Mean Veer (Using Monin-Obukhov approach)
    - Reynolds Stress ``u'w'`` (using Lewkowicz relations)
    - Reynolds Stress ``u'u'``, ``u'v'``, ``u'w'``, ``v'v'``, ``v'w'``, ``w'w'`` (using the :doc:`adem`)
    - Spectra ``Sij`` (using Kaimal, von Karman)
    - Spectra ``Sij`` (using :doc:`adem`)
    - Integral turbulent intensity and lengthscale ``I``, ``l``
#. **Best fit parameter sets**
    - To describe the above profiles analytically (given measured velocity data from an instrument).

In future, generation of artificial flow fields for simulation purposes might be considered. This would overlap with - or replace - utilities like TurbSim to produce, for example:
    - ``.wnd`` fields input to BEM or FVM models like Bladed, FAST and TurbineGRID
    - Inlet boundary conditions for DES or LES codes


Uses
=====

At `Octue <http://www.octue.com>`_, **es-flow** is used to:

  * Provide a basis for developing and validating new processes, like the :doc:`adem`, for characterising Atmospheric Boundary Layers.
  * Process LiDAR and ADCP datasets from raw instrument files.
  * Apply windowed analyses to time series data, for flow characterisation.
  * Generate load cases for FAST, Bladed and our own TurbineGRID aerodynamic wind turbine analysis tools.

We'd like to hear about your use case. Please get in touch!

We use the `GitHub Issue Tracker <https://github.com/octue/es-flow>`_ to manage bug reports and feature requests.



.. toctree::
   :maxdepth: 1
   :hidden:

   self
   models
   examples
   file_formats
   installation
   license
   version_history
   bibliography
   library_api/library_root



.. rubric:: Footnotes

.. [#f1] Not all of this functionality is implemented yet!
