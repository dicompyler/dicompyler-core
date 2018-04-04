#!/usr/bin/env python
# -*- coding: utf-8 -*-
# config.py
"""Configuration for dicompyler-core."""
# Copyright (c) 2016 Aditya Panchal
# This file is part of dicompyler-core, released under a BSD license.
#    See the file license.txt included with this distribution, also
#    available at https://github.com/dicompyler/dicompyler-core/

pil_available = True
try:
    from PIL import Image
except:
    pil_available = False
shapely_available = True
try:
    from shapely.geometry import Polygon
except:
    shapely_available = False
skimage_available = True
try:
    from skimage.transform import rescale
except ImportError:
    skimage_available = False
