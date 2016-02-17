dicompyler-core
===============

|pypi| |travis-ci| |coveralls| |Documentation Status|

Core functionality of `dicompyler <http://www.dicompyler.com>`__. This
package includes:

-  ``dicomparser``: class that parses DICOM objects in an easy-to-use
   manner
-  ``dvhcalc``: independent dose volume histogram (DVH) calculation if dose grid and structure data is present

Other information
-----------------

-  Free software: `BSD license <https://github.com/dicompyler/dicompyler-core/blob/master/LICENSE>`__
-  Documentation: `Read the
   docs <https://dicompyler-core.readthedocs.org>`__
-  Tested on Python 2.7/3.3+

Dependencies
------------

-  `numpy <http://www.numpy.org>`__ 1.2 or higher
-  `pydicom <http://www.pydicom.org>`__ 0.9.9 or higher

   -  pydicom 1.0 is preferred and can be installed via pip using: ``pip install https://github.com/darcymason/pydicom/archive/master.zip``

-  `six <https://pythonhosted.org/six/>`__ 1.5 or higher
-  Optional:

   -  `Pillow <http://python-pillow.org/>`__ (for image display)

Basic Usage
------------

.. code-block:: python

	from dicompylercore import dicomparser
	dp = dicomparser.DicomParser(filename="rtss.dcm")

	# i.e. Get a dict of structure information
	structures = dp.GetStructures()

	>>> structures[5]
	{'color': array([255, 128, 0]), 'type': 'ORGAN', 'id': 5, 'empty': False, 'name': 'Heart'}

Credits
-------

This package was created with
`Cookiecutter <https://github.com/audreyr/cookiecutter>`__ and the
`audreyr/cookiecutter-pypackage <https://github.com/audreyr/cookiecutter-pypackage>`__ project template.

.. |pypi| image:: https://img.shields.io/pypi/v/dicompyler-core.svg
   :target: https://pypi.python.org/pypi/dicompyler-core
.. |travis-ci| image:: https://img.shields.io/travis/dicompyler/dicompyler-core.svg
   :target: https://travis-ci.org/dicompyler/dicompyler-core
.. |coveralls| image:: https://coveralls.io/repos/github/dicompyler/dicompyler-core/badge.svg?branch=master
   :target: https://coveralls.io/github/dicompyler/dicompyler-core?branch=master
.. |Documentation Status| image:: https://readthedocs.org/projects/dicompyler-core/badge/?version=latest
   :target: https://readthedocs.org/projects/dicompyler-core/?badge=latest
