.. _chapter-version-history:

===============
Version History
===============

Origins
=======

EnvironmentSTUDIO flow began as an internal tool at Octue, for characterising turbulence at tidal power sites.
It was originated in MATLAB (then called the 'Flow Characterisation Suite'), and grew incrementally.

Unfortunately, MATLAB is extremely unwieldy if deployed in production, requiring either (complicated and badly limited)
cross compilation, expensive cloud sever licenses or build and deployment with MATLAB's 3Gb+ runtime library
(prohibitively big for practical use in cloud services).

As the move was made toward object orientation, and usage increased toward requiring a production version that could be
deployed without requiring expensive MATLAB licenses, EnvironmentSTUDIO flow project was started in C++.

Gradually, capability is being ported from MATLAB (the original repo having been subsumed within this project and now
being progressively deprecated).


0.1.0
======

This version bump was funded by AURA via the University of Hull. The objective of the work was to validate capabilities
using measured data, and the library was refactored to allow this work to happen in collaboration with Univ. Hull and
the Offshore Renewable Energy Catapult.

New Features
------------
#. Scope of library reduced to preclude specific instrument data readers, and focus on environmental characterisation.
#. Library API consolidated.
#. ADEM functionality ported from legacy MATLAB and tested to work with C++.
#. Added plotting routines (using cpplot).
#. Documentation and build systems implemented and settled.
#. Google Test used to create test harness for the majority of the library.

Backward Incompatible API Changes
---------------------------------
#. Entire API altered to reflect change in scope; will be much more stable going forward.

Bug Fixes & Minor Changes
-------------------------
n/a


0.0.1
======

Initial library release - development version. Highly unstable!

The work to develop this stable release and CI flow undertaken to support research funded by AURA via the University of Hull (see upcoming version 0.1.0).

New Features
------------
#. C++ based `flow` library, with `cmake` build system
#. Large chunks of WIP for 0.1.0 release, included here because the work in progress was used to generate the doc builds.
#. Octue's legacy `es-flow` MATLAB library (to be deprecated)
#. Outline documentation using ReadTheDocs
#. Autodoc build system (with automated build for releases in travis-ci)

Backward Incompatible API Changes
---------------------------------
#. n/a (Initial release)

Bug Fixes & Minor Changes
-------------------------
#. n/a (Initial Release)

