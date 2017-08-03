#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __init__.py
"""setup.py for dicompyler-core."""
# Copyright (c) 2009-2016 Aditya Panchal
# This file is part of dicompyler-core, relased under a BSD license.
#    See the file license.txt included with this distribution, also
#    available at https://github.com/dicompyler/dicompyler-core/

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='dicompyler-core',
    version='0.5.3',
    description="A library of core radiation therapy modules for DICOM / " +
                "DICOM RT used by dicompyler",
    long_description=readme + '\n\n' + history,
    author="Aditya Panchal",
    author_email='apanchal@bastula.org',
    url='https://github.com/dicompyler/dicompyler-core',
    packages=[
        'dicompylercore',
    ],
    package_dir={'dicompylercore':
                 'dicompylercore'},
    include_package_data=True,
    install_requires=[
        "numpy>=1.2",
        "six>=1.5",
        "pydicom>=0.9.9",
        "matplotlib>=1.3.0, <=1.5.1"
    ],
    extras_require={
        'image': ["pillow>=1.0"]
    },
    dependency_links=[
        "git+https://github.com/darcymason/pydicom.git#egg=pydicom-1.0.0"
    ],
    license="BSD License",
    zip_safe=False,
    keywords=[
        'dicompyler-core',
        'dicompylercore',
        'dicompyler'
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: Scientific/Engineering :: Physics'
    ],
    test_suite='tests',
    tests_require=test_requirements
)
