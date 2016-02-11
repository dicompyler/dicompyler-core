dicompyler-core
===============

[![pypi](https://img.shields.io/pypi/v/dicompyler-core.svg)](https://pypi.python.org/pypi/dicompyler-core)

[![travis-ci](https://img.shields.io/travis/dicompyler/dicompyler-core.svg)](https://travis-ci.org/dicompyler/dicompyler-core)

[![Documentation Status](https://readthedocs.org/projects/dicompyler-core/badge/?version=latest)](https://readthedocs.org/projects/dicompyler-core/?badge=latest)

Core functionality of [dicompyler](http://www.dicompyler.com). This package includes:

* `dicomparser`: class that parses DICOM objects in an easy-to-use manner
* DVH calculation: independent dose volume histogram (DVH) calculation if dose grid and structure data is present


Other information
-----------------
* Free software: [BSD license](LICENSE)
* Documentation: [Read the docs](https://dicompyler-core.readthedocs.org)
* Tested on Python 2.7/3.4+

Dependencies
------------

* [numpy](http://www.numpy.org) 1.2 or higher
* [pydicom](http://www.pydicom.org) 1.0 or higher
* [six](https://pythonhosted.org/six/) 1.5 or higher
* Optional:
  * [Pillow](http://python-pillow.org/) (for image display)

Credits
-------

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [`audreyr/cookiecutter-pypackage`](https://github.com/audreyr/cookiecutter-pypackage) project template.
