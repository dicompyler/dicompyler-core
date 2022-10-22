dicompyler-core
===============

|Binder| |pypi| |Python Version| |GH Actions| |Documentation Status| |coveralls| |Codacy| |Codecov| |LGTM| |Total Lines| |Code Size| |Zenodo|

A library of core radiation therapy modules for DICOM / DICOM RT used by `dicompyler <http://www.dicompyler.com>`__. This
package includes:

-  ``dicomparser``: parse DICOM objects in an easy-to-use manner
-  ``dvh``: Pythonic access to dose volume histogram (DVH) data
-  ``dvhcalc``: Independent DVH calculation using DICOM RT Dose & RT Structure Set
-  ``dose``: Pythonic access to RT Dose data including dose summation

Other information
-----------------

-  Free software: `BSD license <https://github.com/dicompyler/dicompyler-core/blob/master/LICENSE>`__
-  Documentation: `Read the docs <https://dicompyler-core.readthedocs.io>`__
-  Tested on Python 3.7+

Dependencies
------------

-  `numpy <http://www.numpy.org>`__ 1.2 or higher
-  `pydicom <https://pydicom.github.io>`__ 0.9.9 or higher (pydicom 1.0 compatible)
-  `matplotlib <http://matplotlib.org>`__ 1.3.0 or higher (for DVH calculation)
-  `six <https://pythonhosted.org/six/>`__ 1.5 or higher
-  Optional:

   -  `Pillow <https://pillow.readthedocs.io>`__ (for image display)
   -  `Shapely <https://github.com/Toblerity/Shapely>`__ (for structure volume calculation)
   -  `scikit-image <http://scikit-image.org/>`__ (for DVH interpolation)
   -  `scipy <https://scipy.org/>`__ (for dose grid summation using interpolation)

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

Citing dicompyler-core
----------------------
A DOI for dicompyler-core with various citation styles can be found at Zenodo: |Zenodo|


Credits
-------

This package was created with
`Cookiecutter <https://github.com/audreyr/cookiecutter>`__ and the
`audreyr/cookiecutter-pypackage <https://github.com/audreyr/cookiecutter-pypackage>`__ project template.

.. |Binder| image:: http://mybinder.org/badge.svg
   :target: http://mybinder.org/repo/bastula/dicom-notebooks
.. |pypi| image:: https://img.shields.io/pypi/v/dicompyler-core.svg
   :target: https://pypi.python.org/pypi/dicompyler-core
.. |Python Version| image:: https://img.shields.io/badge/python-3.7+-blue.svg
   :target: https://pypi.python.org/pypi/dicompyler-core
.. |GH Actions| image:: https://github.com/dicompyler/dicompyler-core/actions/workflows/build.yml/badge.svg
   :target: https://github.com/dicompyler/dicompyler-core/actions
.. |Documentation Status| image:: https://readthedocs.org/projects/dicompyler-core/badge/?version=latest
   :target: https://dicompyler-core.readthedocs.io/en/latest/
.. |coveralls| image:: https://coveralls.io/repos/github/dicompyler/dicompyler-core/badge.svg?branch=master
   :target: https://coveralls.io/github/dicompyler/dicompyler-core?branch=master
.. |Codacy| image:: https://api.codacy.com/project/badge/Grade/27ebb3802baf4d96b0783a2ae5904264
   :target: https://app.codacy.com/gh/dicompyler/dicompyler-core/dashboard
.. |Codecov| image:: https://codecov.io/gh/dicompyler/dicompyler-core/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/dicompyler/dicompyler-core
.. |LGTM| image:: https://img.shields.io/lgtm/alerts/g/dicompyler/dicompyler-core.svg?logo=lgtm&logoWidth=18
   :target: https://lgtm.com/projects/g/dicompyler/dicompyler-core/alerts/
.. |Total Lines| image:: https://img.shields.io/tokei/lines/github/dicompyler/dicompyler-core
   :target: https://img.shields.io/tokei/lines/github/dicompyler/dicompyler-core
.. |Code Size| image:: https://img.shields.io/github/languages/code-size/dicompyler/dicompyler-core
   :target: https://img.shields.io/github/languages/code-size/dicompyler/dicompyler-core
.. |Zenodo| image:: https://zenodo.org/badge/51550203.svg
   :target: https://zenodo.org/badge/latestdoi/51550203
