#!/usr/bin/env python
# -*- coding: utf-8 -*-
# dvhcalc.py
"""Calculate dose volume histogram (DVH) from DICOM RT Structure/Dose data."""
# Copyright (c) 2011-2016 Aditya Panchal
# Copyright (c) 2010 Roy Keyes
# This file is part of dicompyler-core, released under a BSD license.
#    See the file license.txt included with this distribution, also
#    available at https://github.com/dicompyler/dicompyler-core/

from __future__ import division
import numpy as np
import numpy.ma as ma
import matplotlib.path
from dicompylercore import dvh
import collections
from six import iteritems
import logging
logger = logging.getLogger('dicompylercore.dvhcalc')


def get_dvh(structure,
            dose,
            roi,
            limit=None,
            calculate_full_volume=True,
            callback=None):
    """Calculate a cumulative DVH in Gy from a DICOM RT Structure Set & Dose.

    Parameters
    ----------
    structure : pydicom Dataset
        DICOM RT Structure Set used to determine the structure data.
    dose : pydicom Dataset
        DICOM RT Dose used to determine the dose grid.
    roi : int
        The ROI number used to uniquely identify the structure in the structure
        set.
    limit : int, optional
        Dose limit in cGy as a maximum bin for the histogram.
    calculate_full_volume : bool, optional
        Calculate the full structure volume including contours outside of the
        dose grid.
    callback : function, optional
        A function that will be called at every iteration of the calculation.
    """
    from dicompylercore import dicomparser
    rtss = dicomparser.DicomParser(structure)
    rtdose = dicomparser.DicomParser(dose)
    structures = rtss.GetStructures()
    s = structures[roi]
    s['planes'] = rtss.GetStructureCoordinates(roi)
    s['thickness'] = rtss.CalculatePlaneThickness(s['planes'])

    calcdvh = calculate_dvh(s, rtdose, limit, calculate_full_volume, callback)
    return dvh.DVH(counts=calcdvh.histogram,
                   bins=(np.arange(0, 2) if (calcdvh.histogram.size == 1) else
                         np.arange(0, calcdvh.histogram.size + 1) / 100),
                   dvh_type='differential',
                   dose_units='gy',
                   notes=calcdvh.notes,
                   name=s['name']).cumulative


def calculate_dvh(structure,
                  dose,
                  limit=None,
                  calculate_full_volume=True,
                  callback=None):
    """Calculate the differential DVH for the given structure and dose grid.

    Parameters
    ----------
    structure : dict
        A structure (ROI) from an RT Structure Set parsed using DicomParser
    dose : DicomParser
        A DicomParser instance of an RT Dose
    limit : int, optional
        Dose limit in cGy as a maximum bin for the histogram.
    calculate_full_volume : bool, optional
        Calculate the full structure volume including contours outside of the
        dose grid.
    callback : function, optional
        A function that will be called at every iteration of the calculation.
    """
    planes = collections.OrderedDict(sorted(iteritems(structure['planes'])))
    calcdvh = collections.namedtuple('DVH', ['notes', 'histogram'])
    logger.debug("Calculating DVH of %s %s", structure['id'],
                 structure['name'])

    # Create an empty array of bins to store the histogram in cGy
    # only if the structure has contour data or the dose grid exists
    if ((len(planes)) and ("PixelData" in dose.ds)):

        # Get the dose and image data information
        dd = dose.GetDoseData()
        id = dose.GetImageData()

        # Generate a 2d mesh grid to create a polygon mask in dose coordinates
        # Code taken from Stack Overflow Answer from Joe Kington:
        # https://stackoverflow.com/q/3654289/74123
        # Create vertex coordinates for each grid cell
        x, y = np.meshgrid(np.array(dd['lut'][0]), np.array(dd['lut'][1]))
        x, y = x.flatten(), y.flatten()
        dosegridpoints = np.vstack((x, y)).T

        maxdose = int(dd['dosemax'] * dd['dosegridscaling'] * 100)
        # Remove values above the limit (cGy) if specified
        if isinstance(limit, int):
            if (limit < maxdose):
                maxdose = limit
        hist = np.zeros(maxdose)
    else:
        return calcdvh('Empty DVH', np.array([0]))

    n = 0
    notes = None
    planedata = {}
    # Iterate over each plane in the structure
    for z, plane in iteritems(planes):
        # Get the dose plane for the current structure plane
        doseplane = dose.GetDoseGrid(z)
        if doseplane.size:
            planedata[z] = calculate_plane_histogram(plane, doseplane,
                                                     dosegridpoints, maxdose,
                                                     dd, id, structure, hist)
        else:
            # If the dose plane is not found, still perform the calculation
            # but only use it to calculate the volume for the slice
            if not calculate_full_volume:
                logger.warning('Dose plane not found for %s. Contours' +
                               ' not used for volume calculation.', z)
                notes = 'Dose grid does not encompass every contour.' + \
                    ' Volume calculated within dose grid.'
            else:
                origin_z = id['position'][2]
                logger.warning('Dose plane not found for %s.' +
                               ' Using %s to calculate contour volume.',
                               z, origin_z)
                _, vol = calculate_plane_histogram(
                    plane, dose.GetDoseGrid(origin_z), dosegridpoints, maxdose,
                    dd, id, structure, hist)
                planedata[z] = (np.array([0]), vol)
                notes = 'Dose grid does not encompass every contour.' + \
                    ' Volume calculated for all contours.'
        n += 1
        if callback:
            callback(n, len(planes))
    # Volume units are given in cm^3
    volume = sum([p[1] for p in planedata.values()]) / 1000
    # Rescale the histogram to reflect the total volume
    hist = sum([p[0] for p in planedata.values()])
    if hist.max() > 0:
        hist = hist * volume / sum(hist)
    else:
        return calcdvh('Empty DVH', np.array([0]))
    # Remove the bins above the max dose for the structure
    hist = np.trim_zeros(hist, trim='b')

    return calcdvh(notes, hist)


def calculate_plane_histogram(plane, doseplane, dosegridpoints, maxdose, dd,
                              id, structure, hist):
    """Calculate the DVH for the given plane in the structure."""
    contours = [[x[0:2] for x in c['data']] for c in plane]

    # Create a zero valued bool grid
    grid = np.zeros((dd['rows'], dd['columns']), dtype=np.uint8)

    # Calculate the dose plane mask for each contour in the plane
    # and boolean xor to remove holes
    for i, contour in enumerate(contours):
        m = get_contour_mask(dd, id, dosegridpoints, contour)
        grid = np.logical_xor(m.astype(np.uint8), grid).astype(np.bool)

    hist, vol = calculate_contour_dvh(grid, doseplane, maxdose, dd, id,
                                      structure)
    return (hist, vol)


def get_contour_mask(dd, id, dosegridpoints, contour):
    """Get the mask for the contour with respect to the dose plane."""
    doselut = dd['lut']

    c = matplotlib.path.Path(list(contour))
    grid = c.contains_points(dosegridpoints)
    grid = grid.reshape((len(doselut[1]), len(doselut[0])))

    return grid


def calculate_contour_dvh(mask, doseplane, maxdose, dd, id, structure):
    """Calculate the differential DVH for the given contour and dose plane."""
    # Multiply the structure mask by the dose plane to get the dose mask
    mask = ma.array(doseplane * dd['dosegridscaling'] * 100, mask=~mask)
    # Calculate the differential dvh
    hist, edges = np.histogram(mask.compressed(),
                               bins=maxdose,
                               range=(0, maxdose))

    # Calculate the volume for the contour for the given dose plane
    vol = sum(hist) * ((id['pixelspacing'][0]) *
                       (id['pixelspacing'][1]) *
                       (structure['thickness']))
    return hist, vol
