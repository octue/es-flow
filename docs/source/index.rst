=================
EnvironmentSTUDIO
=================

EnvironmentSTUDIO flow (**es-flow** for short) is a collection of libraries and executables
for analysing and modeling both atmospheric and marine boundary layers.
es-flow covers [#f1]_ the following steps:

a. Read-in and clean up raw instrument data from a variety of instruments including LiDAR, SODAR, Met Masts, ADCPs, Vectrinos and Microstructure Profilers. For a full list of supported instruments and operational modes see For more, see :doc:`supported_instruments`.

b. Save and load cleaned instrument data in manageable, structured file formats based on OAS standards. For more, see :doc:`file_formats`.

c. Apply conventional flow analyses for wind and tidal site characterisation using instrument data.

d. Apply spectral and coherent-structural turbulence analyses for wind and tidal site characterisation.

e. Save analysis results files in OAS standard forms (again, see :doc:`file_formats`).

f. Generate artificial flow fields for simulation purposes (fields for FVM/panel/BEM models, inlets to DES, etc.).

g. Provide an API for accessing artificial flow field properties (e.g. velocity fields) and/or saving the generated fields to OAS standard or industry standard files.

We use the `GitHub Issue Tracker <https://github.com/oceanarraysystems/es-flow>`_
to manage bug reports and feature requests.


.. toctree::
   :maxdepth: 1
   :hidden:

   installation
   api/library_root
   license
   users
   supported_instruments
   file_formats
   examples
   version_history
   bibliography



.. rubric:: Footnotes

.. [#f1] Or **will** cover. Not all of this functionality is implemented yet!
