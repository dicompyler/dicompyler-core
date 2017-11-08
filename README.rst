dicompyler-core
===============

|Binder| |pypi| |travis-ci| |circleci| |coveralls| |Codacy| |Codecov| |Documentation Status|

Core functionality of `dicompyler <http://www.dicompyler.com>`__. This
package includes:

-  ``dicomparser``: parse DICOM objects in an easy-to-use manner
-  ``dvh``: Pythonic access to dose volume histogram (DVH) data
-  ``dvhcalc``: independent dose volume histogram (DVH) calculation if dose grid and structure data is present

Other information
-----------------

-  Free software: `BSD license <https://github.com/dicompyler/dicompyler-core/blob/master/LICENSE>`__
-  Documentation: `Read the
   docs <https://dicompyler-core.readthedocs.org>`__
-  Tested on Python 2.7/3.4+

Dependencies
------------

-  `numpy <http://www.numpy.org>`__ 1.2 or higher
-  `pydicom <http://www.pydicom.org>`__ 0.9.9 or higher

   -  pydicom 1.0 is preferred and can be installed via pip using: ``pip install https://github.com/darcymason/pydicom/archive/master.zip --upgrade``

-  `matplotlib <http://matplotlib.org>`__ 1.3.0 or higher (for DVH calculation)
-  `six <https://pythonhosted.org/six/>`__ 1.5 or higher
-  Optional:

   -  `Pillow <http://python-pillow.org/>`__ (for image display)
   -  `Shapely <https://github.com/Toblerity/Shapely>`__ (for structure volume calculation)

Basic Usage
------------

.. code-block:: python

    from dicompylercore import dicomparser, dvh, dvhcalc
    dp = dicomparser.DicomParser("rtss.dcm")

    # i.e. Get a dict of structure information
    structures = dp.GetStructures()

    >>> structures[5]
    {'color': array([255, 128, 0]), 'type': 'ORGAN', 'id': 5, 'empty': False, 'name': 'Heart'}

    # Access DVH data
    rtdose = dicomparser.DicomParser("rtdose.dcm")
    heartdvh = dvh.DVH.from_dicom_dvh(rtdose.ds, 5)

    >>> heartdvh.describe()
    Structure: Heart
    DVH Type:  cumulative, abs dose: Gy, abs volume: cm3
    Volume:    437.46 cm3
    Max Dose:  3.10 Gy
    Min Dose:  0.02 Gy
    Mean Dose: 0.64 Gy
    D100:      0.00 Gy
    D98:       0.03 Gy
    D95:       0.03 Gy
    D2cc:      2.93 Gy

    # Calculate a DVH from DICOM RT data
    calcdvh = dvhcalc.get_dvh("rtss.dcm", "rtdose.dcm", 5)

    >>> calcdvh.max, calcdvh.min, calcdvh.D2cc
    (3.0899999999999999, 0.029999999999999999, dvh.DVHValue(2.96, 'Gy'))

Advanced Usage and Examples can be found in Binder: |Binder|

Credits
-------

This package was created with
`Cookiecutter <https://github.com/audreyr/cookiecutter>`__ and the
`audreyr/cookiecutter-pypackage <https://github.com/audreyr/cookiecutter-pypackage>`__ project template.

.. |Binder| image:: http://mybinder.org/badge.svg
   :target: http://mybinder.org/repo/bastula/dicom-notebooks
.. |pypi| image:: https://img.shields.io/pypi/v/dicompyler-core.svg
   :target: https://pypi.python.org/pypi/dicompyler-core
.. |travis-ci| image:: https://img.shields.io/travis/dicompyler/dicompyler-core.svg
   :target: https://travis-ci.org/dicompyler/dicompyler-core
.. |circleci| image:: https://circleci.com/gh/dicompyler/dicompyler-core.svg?style=svg
   :target: https://circleci.com/gh/dicompyler/dicompyler-core
.. |coveralls| image:: https://coveralls.io/repos/github/dicompyler/dicompyler-core/badge.svg?branch=master
   :target: https://coveralls.io/github/dicompyler/dicompyler-core?branch=master
.. |Codacy| image:: https://api.codacy.com/project/badge/Grade/d8e948ed96914dc19293e34060847d3f
   :target: https://www.codacy.com/app/bastula/dicompyler-core?utm_campaign=Badge_Coverage
.. |Codecov| image:: https://codecov.io/gh/dicompyler/dicompyler-core/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/dicompyler/dicompyler-core
.. |Documentation Status| image:: https://readthedocs.org/projects/dicompyler-core/badge/?version=latest
   :target: https://dicompyler-core.readthedocs.io/en/latest/
