=====
Usage
=====

To use dicompyler-core in a project:

.. code-block:: python

    from dicompylercore import dicomparser
    dp = dicomparser.DicomParser(filename="rtss.dcm")

    # i.e. Get a dict of structure information
    structures = dp.GetStructures()

    >>> structures[5]
    {'color': array([255, 128, 0]), 'type': 'ORGAN', 'id': 5, 'empty': False, 'name': 'Heart'}