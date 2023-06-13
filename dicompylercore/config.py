#!/usr/bin/env python
# -*- coding: utf-8 -*-
# config.py
"""Configuration for dicompyler-core."""
# Copyright (c) 2016-2018 Aditya Panchal
# This file is part of dicompyler-core, released under a BSD license.
#    See the file license.txt included with this distribution, also
#    available at https://github.com/dicompyler/dicompyler-core/

mpl_available = True
pil_available = True
shapely_available = True
skimage_available = True
scipy_available = True

import importlib

mpl_available = importlib.util.find_spec("matplotlib") is not None
pil_available = importlib.util.find_spec("PIL") is not None
shapely_available = importlib.util.find_spec("shapely") is not None
skimage_available = importlib.util.find_spec("skimage") is not None
scipy_available = importlib.util.find_spec("scipy") is not None


# DICOM UID prefix
dicompyler_uid_prefix = "1.2.826.0.1.3680043.8.1070."
dicompyler_uid_prefix_image = dicompyler_uid_prefix + "1."
dicompyler_uid_prefix_rtstruct = dicompyler_uid_prefix + "2."
dicompyler_uid_prefix_rtplan = dicompyler_uid_prefix + "3."
dicompyler_uid_prefix_rtdose = dicompyler_uid_prefix + "4."
