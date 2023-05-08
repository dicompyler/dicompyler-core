=======
History
=======

0.5.7 (unreleased)
------------------

0.5.6 (2023-05-08)
------------------

.. warning:: This version will be the last version to support Python 2.x and support will be dropped in version 0.5.7.
- Dropped support for Python 3.5 & 3.6 and added support for Python 3.9 & 3.10.
- Made changes to codebase to support recent versions of numpy, Shapely and scikit-image dependencies.
- Added ``dose`` module with ``DVH`` class for Pythonic access to RT Dose. `(#164) <https://github.com/dicompyler/dicompyler-core/pull/164>`__ `[Dan Cutright] <https://github.com/cutright>`__
- Added decubitus orientation and related changes. `(#285) <https://github.com/dicompyler/dicompyler-core/pull/285>`__ `[Darcy Mason] <https://github.com/darcymason>`__
- Fix a bug if Pixel Data attribute was set for non image based SOP Classes (i.e. RT Structure Set). `(#214) <https://github.com/dicompyler/dicompyler-core/pull/214>`__ `[Dan Cutright] <https://github.com/cutright>`__

dvhcalc
~~~~~~~
- Implement interpolation for non square pixels in DVH calculation. `(#124) <https://github.com/dicompyler/dicompyler-core/pull/124>`__
- Fix a bug where the DVHDoseScaling attribute was not applied properly to RT Dose DVHs. `(#301) <https://github.com/dicompyler/dicompyler-core/pull/301>`__ `[Christian Velten] <https://github.com/cvelten>`__
- Fix a bug where floating point pixel spacing wasn't rounded in DVH calculations. `(#318) <https://github.com/dicompyler/dicompyler-core/pull/318>`__ `[Samuel Ouellet] <https://github.com/smichi23>`__

dose
~~~~~~~
- Added RT Dose grid summmation with interpolation (from DVHA). `(#164) <https://github.com/dicompyler/dicompyler-core/pull/164>`__ `[Dan Cutright] <https://github.com/cutright>`__

dicomparser
~~~~~~~~~~~
- Initial implementation of memory mapped access to pixel data. `(#131) <https://github.com/dicompyler/dicompyler-core/pull/131>`__
- Ensure that all files read have a valid File Meta header.

0.5.5 (2019-05-31)
------------------

dvhcalc
~~~~~~~
- Refactored bounding & resampling set up code to only execute
  if conditions are met.
- Fix a bug where the resampled LUT was not calculated
  correctly for DVH interpolation.

dvh
~~~
- Differential DVH calculation modified. `(#60) <https://github.com/dicompyler/dicompyler-core/pull/60>`__ `[Hideki Nakamoto] <https://github.com/inamoto85>`__
- Fix an issue with D100 not returning 0 Gy. `(#74) <https://github.com/dicompyler/dicompyler-core/pull/74>`__ `[Gabriel Couture] <https://github.com/gacou54>`__
- Preserve global maximum dose. `(#106) <https://github.com/dicompyler/dicompyler-core/pull/106>`__ `[Akihisa Wakita] <https://github.com/wkt84>`__

dicomparser
~~~~~~~~~~~
- Remove the test for existence of `ContourImageSequence` in
  `GetStructureCoordinates`. `(#81) <https://github.com/dicompyler/dicompyler-core/pull/81>`__ `[Gabriel Couture] <https://github.com/gacou54>`__
- Utilize integer division when generating a background for
  an image.
- Return a string for the patient's name as `PersonName3`
  cannot be serialized.
- Fix a bug to return a referenced FoR if the
  FrameOfReference is blank.
- Fix a bug in `GetPlan` where the wrong object names were
  used. `(#43) <https://github.com/dicompyler/dicompyler-core/pull/43>`__ `[gertsikkema] <https://github.com/gertsikkema>`__
- Ensure that Rx Dose from RT Plan is rounded instead of
  truncated.
- Account for holes and bifurcated structures for structure
  volume calculation.
- Implement structure volume calculation using Shapely.
- Fix window calculation if not present in header.
- Add checks in max, mean, min and dose_constraint for case where counts array is empty or all 0's. `(#96) <https://github.com/dicompyler/dicompyler-core/pull/96>`__ `[Nicolas Galler] <https://github.com/nicocrm>`__


0.5.4 (2018-04-02)
------------------

dvhcalc
~~~~~~~
- Implemented DVH interpolation. (#39)
- Implemented optional user-specified structure thickness
  for DVH calculation.


dvh
~~~
- Fix a bug in absolute_volume if a DVH instance's volume units
  don't use default of Gy.
- Fix a bug in absolute_dose if a DVH instance's dose units don't
  use default of Gy. (#19)
- Support decimal values for volume constraints (i.e.V71.6).
- Support decimal values for dose constraints (i.e. D0.03cc).

dicomparser
~~~~~~~~~~~
- Ensure that Rx Dose from RT Plan is rounded instead of
  truncated.
- Account for holes and bifurcated structures for structure
  volume calculation.
- Implement structure volume calculation using Shapely. (#28)


0.5.3 (2017-08-03)
------------------
* Added support for plotting structure colors.
* Support Python 2 unicode filenames in dicomparser.
* Support DVH calculation of structures partially covered by the dose grid.


0.5.2 (2016-07-25)
------------------

* Added ``DVH`` class for Pythonic access to dose volume histogram data.
* Refactored and added unit tests for dvhcalc.
* Added examples and usage for ``dvh`` and ``dvhcalc`` modules.
* Jupyter notebook of examples can be found in Binder: |dicom-notebooks|


0.5.1 (2016-02-17)
------------------

* Added support for pydicom 0.9.9 so releases from PyPI can be built.


0.5.0 (2016-02-11)
------------------

* First release on PyPI.

.. |dicom-notebooks| image:: http://mybinder.org/badge.svg
   :target: http://mybinder.org/repo/bastula/dicom-notebooks
