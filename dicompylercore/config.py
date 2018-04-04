#!/usr/bin/env python
# -*- coding: utf-8 -*-
# config.py
"""Configuration for dicompyler-core."""
# Copyright (c) 2016-2018 Aditya Panchal
# This file is part of dicompyler-core, released under a BSD license.
#    See the file license.txt included with this distribution, also
#    available at https://github.com/dicompyler/dicompyler-core/

from six import PY2

pil_available = True
shapely_available = True
skimage_available = True

if PY2:
    import imp
    try:
        imp.find_module('PIL')
    except ImportError:
        pil_available = False

    try:
        imp.find_module('shapely')
    except ImportError:
        shapely_available = False

    try:
        imp.find_module('skimage')
    except ImportError:
        skimage_available = False
else:
    import importlib
    pil_available = importlib.util.find_spec('PIL') is not None
    shapely_available = importlib.util.find_spec('shapely') is not None
    skimage_available = importlib.util.find_spec('skimage') is not None
