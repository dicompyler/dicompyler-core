=====
Usage
=====

To use dicompyler-core in a project:

DICOM data can be easily accessed using convenience functions using the :mod:`dicompylercore.dicomparser.DicomParser` class:

.. code-block:: python

    from dicompylercore import dicomparser, dvh, dvhcalc
    dp = dicomparser.DicomParser("rtss.dcm")

    # i.e. Get a dict of structure information
    structures = dp.GetStructures()

    >>> structures[5]
    {'color': array([255, 128, 0]), 'type': 'ORGAN', 'id': 5, 'empty': False, 'name': 'Heart'}

Dose volume histogram (DVH) data can be accessed in a Pythonic manner using the :mod:`dicompylercore.dvh.DVH` class:

.. code-block:: python

    rtdose = dicomparser.DicomParser("rtdose.dcm")
    heartdvh = dvh.DVH.from_dicom_dvh(rtdose.ds, 5)

    >>> heartdvh.describe()
    Structure: Heart
    -----
    DVH Type:  cumulative, abs dose: Gy, abs volume: cm3
    Volume:    437.46 cm3
    Max Dose:  3.10 Gy
    Min Dose:  0.02 Gy
    Mean Dose: 0.64 Gy
    D100:      0.00 Gy
    D98:       0.03 Gy
    D95:       0.03 Gy
    D2cc:      2.93 Gy

    >>> heartdvh.max, heartdvh.min, heartdvh.D2cc
    (3.0999999999999779, 0.02, dvh.DVHValue(2.9299999999999815, 'Gy'))

Dose volume histograms (DVHs) can be independently calculated using the dvh :mod:`dicompylercore.dvhcalc` module:

.. code-block:: python

    # Calculate a DVH from DICOM RT data
    calcdvh = dvhcalc.get_dvh("rtss.dcm", "rtdose.dcm", 5)

    >>> calcdvh.max, calcdvh.min, calcdvh.D2cc
    (3.0899999999999999, 0.029999999999999999, dvh.DVHValue(2.96, 'Gy'))
