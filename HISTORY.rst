=======
History
=======

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
