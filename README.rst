dicompyler-core
===============

|pypi|

|travis-ci|

|Documentation Status|

Core functionality of `dicompyler <http://www.dicompyler.com>`__. This
package includes:

-  ``dicomparser``: class that parses DICOM objects in an easy-to-use
   manner
-  DVH calculation: independent dose volume histogram (DVH) calculation
   if dose grid and structure data is present

Other information
-----------------

-  Free software: `BSD license <LICENSE>`__
-  Documentation: `Read the
   docs <https://dicompyler-core.readthedocs.org>`__
-  Tested on Python 2.7/3.4+

Dependencies
------------

-  `numpy <http://www.numpy.org>`__ 1.2 or higher
-  `pydicom <http://www.pydicom.org>`__ 1.0 or higher
-  `six <https://pythonhosted.org/six/>`__ 1.5 or higher
-  Optional:

  -  `Pillow <http://python-pillow.org/>`__ (for image display)

Credits
-------

This package was created with
`Cookiecutter <https://github.com/audreyr/cookiecutter>`__ and the
`audreyr/cookiecutter-pypackage <https://github.com/audreyr/cookiecutter-pypackage>`__

.. |pypi| image:: https://img.shields.io/pypi/v/dicompyler-core.svg
   :target: https://pypi.python.org/pypi/dicompyler-core
.. |travis-ci| image:: https://img.shields.io/travis/dicompyler/dicompyler-core.svg
   :target: https://travis-ci.org/dicompyler/dicompyler-core
.. |Documentation Status| image:: https://readthedocs.org/projects/dicompyler-core/badge/?version=latest
   :target: https://readthedocs.org/projects/dicompyler-core/?badge=latest
