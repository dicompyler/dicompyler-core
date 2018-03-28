#!/usr/bin/env python
# -*- coding: utf-8 -*-
# dvhcalc.py
"""Calculate dose volume histogram (DVH) from DICOM RT Structure/Dose data."""
# Copyright (c) 2011-2018 Aditya Panchal
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

skimage_available = True
try:
    from skimage.transform import rescale
except ImportError:
    skimage_available = False


def get_dvh(structure,
            dose,
            roi,
            limit=None,
            calculate_full_volume=True,
            use_structure_extents=False,
            interpolation_resolution=None,
            interpolation_segments_between_planes=0,
            thickness=None,
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
    use_structure_extents : bool, optional
        Limit the DVH calculation to the in-plane structure boundaries.
    interpolation_resolution : float, optional
        Resolution in mm to interpolate the structure and dose data to.
    interpolation_segments_between_planes : integer, optional
        Number of segments to interpolate between structure slices.
    thickness : float, optional
        Structure thickness used to calculate volume of a voxel.
    callback : function, optional
        A function that will be called at every iteration of the calculation.
    """
    from dicompylercore import dicomparser
    rtss = dicomparser.DicomParser(structure)
    rtdose = dicomparser.DicomParser(dose)
    structures = rtss.GetStructures()
    s = structures[roi]
    s['planes'] = rtss.GetStructureCoordinates(roi)
    s['thickness'] = thickness if thickness else rtss.CalculatePlaneThickness(
        s['planes'])

    calcdvh = calculate_dvh(s, rtdose, limit, calculate_full_volume,
                            use_structure_extents, interpolation_resolution,
                            interpolation_segments_between_planes,
                            callback)
    return dvh.DVH(counts=calcdvh.histogram,
                   bins=(np.arange(0, 2) if (calcdvh.histogram.size == 1) else
                         np.arange(0, calcdvh.histogram.size + 1) / 100),
                   dvh_type='differential',
                   dose_units='Gy',
                   notes=calcdvh.notes,
                   name=s['name']).cumulative


def calculate_dvh(structure,
                  dose,
                  limit=None,
                  calculate_full_volume=True,
                  use_structure_extents=False,
                  interpolation_resolution=None,
                  interpolation_segments_between_planes=0,
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
    use_structure_extents : bool, optional
        Limit the DVH calculation to the in-plane structure boundaries.
    interpolation_resolution : float, optional
        Resolution in mm to interpolate the structure and dose data to.
    interpolation_segments_between_planes : integer, optional
        Number of segments to interpolate between structure slices.
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

        dgindexextents = []
        if use_structure_extents:
            # Determine structure and respectively dose grid extents
            extents = structure_extents(structure['planes'])
            dgindexextents = dosegrid_extents_indices(extents, dd)
            dgextents = dosegrid_extents_positions(dgindexextents, dd)
            # Determine LUT from extents
            if interpolation_resolution:
                dd['lut'] = lut_from_extents(
                    dgextents, sampling_rate=interpolation_resolution,
                    min_sampling_rate=id['pixelspacing'][0])
                dd['rows'] = dd['lut'][1].shape[0]
                dd['columns'] = dd['lut'][0].shape[0]

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
    # Interpolate between planes in the direction of the structure
    if interpolation_segments_between_planes:
        planes = interpolate_between_planes(
            planes, interpolation_segments_between_planes)
        # Thickness derived from total number of segments relative to original
        structure['thickness'] = structure[
            'thickness'] / (interpolation_segments_between_planes + 1)

    # Iterate over each plane in the structure
    for z, plane in iteritems(planes):
        # Get the dose plane for the current structure plane
        if interpolation_resolution and not skimage_available:
            raise ImportError(
                "scikit-image must be installed to perform DVH interpolation.")
        if interpolation_resolution and skimage_available:
            doseplane = get_interpolated_dose(
                dose, z, interpolation_resolution, dgindexextents)
        else:
            doseplane = dose.GetDoseGrid(z)
        if doseplane.size:
            planedata[z] = calculate_plane_histogram(plane, doseplane,
                                                     dosegridpoints, maxdose,
                                                     dd, id, structure, hist)
            # print(f'Slice: {z}, volume: {planedata[z][1]}')
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
    # print(f'total volume: {volume}')
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

    # def inpolygon(polygon, xp, yp):
    #     return np.array(
    #         [Point(x, y).intersects(polygon) for x, y in zip(xp, yp)],
    #         dtype=np.bool)

    # p = Polygon(contour)
    # x, y = np.meshgrid(np.array(dd['lut'][0]), np.array(dd['lut'][1]))
    # mask = inpolygon(p, x.ravel(), y.ravel())
    # return mask.reshape((len(doselut[1]), len(doselut[0])))

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
    vol = sum(hist) * ((np.mean(np.diff(dd['lut'][0]))) *
                       (np.mean(np.diff(dd['lut'][1]))) *
                       (structure['thickness']))
    return hist, vol


def structure_extents(coords):
    """Determine structure extents in patient coordinates."""
    bounds = []
    for i, z in enumerate(sorted(coords.items())):
        contours = [[x[0:2] for x in c['data']] for c in z[1]]
        for contour in contours:
            x, y = np.array([x[0:1] for x in contour]), np.array(
                [x[1:2] for x in contour])
            bounds.append([np.min(x), np.min(y), np.max(x), np.max(y)])
    extents = np.array(bounds)
    return np.array(
        [np.amin(extents, axis=0)[0:2],
         np.amax(extents, axis=0)[2:4]]).flatten().tolist()


def dosegrid_extents_indices(extents, dd, padding=1):
    """Determine dose grid extents as array indices."""
    dgxmin = np.argmin(np.fabs(dd['lut'][0] - extents[0])) - padding
    if dd['lut'][0][dgxmin] > extents[0]:
        dgxmin -= 1
    dgxmax = np.argmin(np.fabs(dd['lut'][0] - extents[2])) + padding
    dgymin = np.argmin(np.fabs(dd['lut'][1] - extents[1])) - padding
    dgymax = np.argmin(np.fabs(dd['lut'][1] - extents[3])) + padding
    dgxmin = 0 if dgxmin < 0 else dgxmin
    dgymin = 0 if dgymin < 0 else dgymin
    if dgxmax == dd['lut'][0].shape[0]:
        dgxmax = dd['lut'][0].shape[0] - 1
    if dgymax == dd['lut'][1].shape[0]:
        dgymax = dd['lut'][1].shape[0] - 1
    return np.array([dgxmin, dgymin, dgxmax, dgymax])


def dosegrid_extents_positions(extents, dd):
    """Determine dose grid extents in patient coordinate indices."""
    return np.array([dd['lut'][0][extents[0]],
                     dd['lut'][1][extents[1]],
                     dd['lut'][0][extents[2]],
                     dd['lut'][1][extents[3]]])


def lut_from_extents(extents, sampling_rate, min_sampling_rate=3):
    """Determine new patient to pixel LUT from sampling rate."""
    if (min_sampling_rate % sampling_rate != 0.0):
        raise AttributeError(
            "Sampling rate must be a factor of %g/(2^n)," % min_sampling_rate +
            " where n is an integer. Value provided was %g." % sampling_rate)
    xsamples = int(abs((extents[0] - extents[2]) / sampling_rate) + 1)
    ysamples = int(abs((extents[1] - extents[3]) / sampling_rate) + 1)
    x = np.linspace(extents[0], extents[2], xsamples, dtype=np.float)
    y = np.linspace(extents[1], extents[3], ysamples, dtype=np.float)
    return x[:-1], y[:-1]


def get_interpolated_dose(dose, z, resolution, extents=[]):
    """Get interpolated dose for the given z, resolution & array extents."""
    d = dose.GetDoseGrid(z)
    scale = (np.array(dose.ds.PixelSpacing) / resolution).tolist()
    extent_dose = d[extents[1]:extents[3],
                    extents[0]:extents[2]] if len(extents) else d
    interp_dose = rescale(
        extent_dose,
        scale=scale,
        mode='symmetric',
        order=1,
        preserve_range=True)
    return interp_dose


def interpolate_between_planes(planes, n=2):
    """Interpolate n additional structure planes (segments) in between planes.
    """
    keymap = {np.array([k], dtype=np.float32)[0]: k for k in planes.keys()}
    sorted_keys = np.sort(np.array(list(planes.keys()), dtype=np.float32))
    num_new_samples = (len(planes.keys()) * (n + 1)) - n
    newgrid = np.linspace(sorted_keys[0], sorted_keys[-1], num_new_samples)
    new_planes = {}
    # If the plane already exists in the dictionary, use it
    # otherwise use the closest plane
    # TODO: Add actual interpolation of structure data between planes
    for i, plane in enumerate(newgrid):
        new_plane = sorted_keys[np.argmin(np.fabs(sorted_keys - plane))]
        new_planes[plane] = planes[keymap[new_plane]]
    return new_planes
