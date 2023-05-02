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
from dicompylercore.config import skimage_available
import collections
try:
    from collections.abc import Sequence
except ImportError:
    from collections import Sequence
from six import iteritems
import logging
logger = logging.getLogger('dicompylercore.dvhcalc')

if skimage_available:
    from skimage.transform import rescale


def get_dvh(structure,
            dose,
            roi,
            limit=None,
            calculate_full_volume=True,
            use_structure_extents=False,
            interpolation_resolution=None,
            interpolation_segments_between_planes=0,
            thickness=None,
            memmap_rtdose=False,
            callback=None):
    """Calculate a cumulative DVH in Gy from a DICOM RT Structure Set & Dose.

    Parameters
    ----------
    structure : pydicom Dataset or filename
        DICOM RT Structure Set used to determine the structure data.
    dose : pydicom Dataset or filename
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
    interpolation_resolution : tuple or float, optional
        Resolution in mm (row, col) to interpolate structure and dose data to.
        If float is provided, original dose grid pixel spacing must be square.
    interpolation_segments_between_planes : integer, optional
        Number of segments to interpolate between structure slices.
    thickness : float, optional
        Structure thickness used to calculate volume of a voxel.
    memmap_rtdose : bool, optional
        Use memory mapping to access the pixel array of the DICOM RT Dose.
        This reduces memory usage at the expense of increased calculation time.
    callback : function, optional
        A function that will be called at every iteration of the calculation.

    Returns
    -------
    dvh.DVH
        An instance of dvh.DVH in cumulative dose. This can be converted to
        different formats using the attributes and properties of the DVH class.
    """
    from dicompylercore import dicomparser

    rtss = dicomparser.DicomParser(structure)
    rtdose = dicomparser.DicomParser(dose, memmap_pixel_array=memmap_rtdose)
    structures = rtss.GetStructures()
    s = structures[roi]
    s['planes'] = rtss.GetStructureCoordinates(roi)
    s['thickness'] = thickness if thickness else rtss.CalculatePlaneThickness(
        s['planes'])

    calcdvh = _calculate_dvh(s, rtdose, limit, calculate_full_volume,
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


def _calculate_dvh(structure,
                   dose,
                   limit=None,
                   calculate_full_volume=True,
                   use_structure_extents=False,
                   interpolation_resolution=None,
                   interpolation_segments_between_planes=0,
                   callback=None):
    """Calculate a differential DVH for the given structure and dose grid.

    Parameters
    ----------
    structure : dict
        A structure (ROI) from an RT Structure Set parsed using DicomParser.
        The dictionary must include a `thickness` key with a thickness `float`.
    dose : DicomParser
        A DicomParser instance of an RT Dose
    limit : int, optional
        Dose limit in cGy as a maximum bin for the histogram.
    calculate_full_volume : bool, optional
        Calculate the full structure volume including contours outside of the
        dose grid.
    use_structure_extents : bool, optional
        Limit the DVH calculation to the in-plane structure boundaries.
    interpolation_resolution : tuple or float, optional
        Resolution in mm (row, col) to interpolate structure and dose data to.
        If float is provided, original dose grid pixel spacing must be square.
    interpolation_segments_between_planes : integer, optional
        Number of segments to interpolate between structure slices.
    callback : function, optional
        A function that will be called at every iteration of the calculation.

    Returns
    -------
    calcdvh: collections.namedtuple
        A named tuple of notes and and histogram data.

    Notes
    -----
    This is an internal function called by `get_dvh` and
    should not be called directly.
    """
    planes = collections.OrderedDict(sorted(iteritems(structure['planes'])))
    calcdvh = collections.namedtuple('DVH', ['notes', 'histogram'])
    logger.debug("Calculating DVH of %s %s", structure['id'],
                 structure['name'])

    # Create an empty array of bins to store the histogram in cGy
    # only if the structure has contour data or the dose grid exists
    if ((len(planes)) and (hasattr(dose, 'pixel_array'))):

        # Get the dose and image data information
        dd = dose.GetDoseData()
        id = dose.GetImageData()

        # Determine structure and respectively dose grid extents
        if interpolation_resolution or use_structure_extents:
            extents = []
            if use_structure_extents:
                extents = structure_extents(structure['planes'])
            dgindexextents = dosegrid_extents_indices(extents, dd)
            dgextents = dosegrid_extents_positions(dgindexextents, dd)
            # Determine LUT from extents
            if use_structure_extents:
                dd['lut'] = \
                    (dd['lut'][0][dgindexextents[0]:dgindexextents[2]],
                     dd['lut'][1][dgindexextents[1]:dgindexextents[3]])
            # If interpolation is enabled, generate new LUT from extents
            if interpolation_resolution:
                dd['lut'] = get_resampled_lut(
                    dgindexextents,
                    dgextents,
                    new_pixel_spacing=interpolation_resolution,
                    min_pixel_spacing=id['pixelspacing'])
            dd['rows'] = dd['lut'][1].shape[0]
            dd['columns'] = dd['lut'][0].shape[0]

        # Generate a 2d mesh grid to create a polygon mask in dose coordinates
        # Code taken from Stack Overflow Answer from Joe Kington:
        # https://stackoverflow.com/q/3654289/74123
        # Create vertex coordinates for each grid cell
        x_index = dd['x_lut_index']
        x, y = np.meshgrid(
            np.array(dd['lut'][x_index]), np.array(dd['lut'][1-x_index])
        )
        x, y = x.flatten(), y.flatten()
        dosegridpoints = np.vstack((x, y)).T

        maxdose = int(dd['dosemax'] * dd['dosegridscaling'] * 100) + 1
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
        if interpolation_resolution or use_structure_extents:
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
                # Create a dummy dose grid with the correct size.
                # calculate_plane_histogram() and its methods provide the
                # volume calc needed, but do so as part of DVH calc,
                # which requires a dose grid
                dummy_dose = dose.GetDoseGrid(origin_z)
                if use_structure_extents:
                    extents = dgindexextents
                    dummy_dose = dummy_dose[
                        extents[1]:extents[3], extents[0]:extents[2]
                    ]
                _, vol = calculate_plane_histogram(
                    plane, dummy_dose, dosegridpoints, maxdose,
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
        grid = np.logical_xor(m.astype(np.uint8), grid).astype(np.bool_)

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
    if dd['x_lut_index'] == 0:  # X values across columns
        grid = grid.reshape((len(doselut[1]), len(doselut[0])))
    else:  # decubitus
        grid = grid.reshape((len(doselut[0]), len(doselut[1]))).T

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
    vol = sum(hist) * (abs(np.mean(np.diff(dd['lut'][0]))) *
                       abs(np.mean(np.diff(dd['lut'][1]))) *
                       (structure['thickness']))
    return hist, vol


def structure_extents(coords):
    """Determine structure extents in patient coordinates.

    Parameters
    ----------
    coords : dict
        Structure coordinates from dicomparser.GetStructureCoordinates.

    Returns
    -------
    list
        Structure extents in patient coordintes: [xmin, ymin, xmax, ymax].
    """
    bounds = []
    for z in sorted(coords.items()):
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
    """Determine dose grid extents from structure extents as array indices.

    Parameters
    ----------
    extents : list
        Structure extents in patient coordinates: [xmin, ymin, xmax, ymax].
        If an empty list, no structure extents will be used in the calculation.
    dd : dict
        Dose data from dicomparser.GetDoseData.
    padding : int, optional
        Pixel padding around the structure extents.

    Returns
    -------
    list
        Dose grid extents in pixel coordinates as array indices:
        [col_min, row_min, col_max, row_max].
    """
    col_lut, row_lut = dd['lut']
    num_cols = len(col_lut)
    num_rows = len(row_lut)
    if not len(extents):
        return [0, 0, num_cols - 1, num_rows - 1]

    if dd['x_lut_index'] == 0:  # X is across rows
        strx_col_min, strx_col_max = extents[0], extents[2]
        strx_row_min, strx_row_max = extents[1], extents[3]
    else:  # decubitus case, X down rows
        strx_col_min, strx_col_max = extents[1], extents[3]
        strx_row_min, strx_row_max = extents[0], extents[2]

    dg_col_min = np.argmin(np.fabs(col_lut - strx_col_min))
    dg_col_max = np.argmin(np.fabs(col_lut - strx_col_max))
    if dg_col_min > dg_col_max:
        dg_col_min, dg_col_max = dg_col_max, dg_col_min
    if col_lut[dg_col_min] > strx_col_min:
        dg_col_min -= 1

    dg_row_min = np.argmin(np.fabs(row_lut - strx_row_min))
    dg_row_max = np.argmin(np.fabs(row_lut - strx_row_max))
    if dg_row_min > dg_row_max:
        dg_row_min, dg_row_max = dg_row_max, dg_row_min

    # Ensure indices within array limits regardless of padding
    dg_col_min = max(0, dg_col_min-padding)
    dg_row_min = max(0, dg_row_min-padding)
    dg_col_max = min(num_cols - 1, dg_col_max+padding)
    dg_row_max = min(num_rows - 1, dg_row_max+padding)

    return [dg_col_min, dg_row_min, dg_col_max, dg_row_max]


def dosegrid_extents_positions(extents, dd):
    """Determine dose grid extents in patient coordinate indices.

    Parameters
    ----------
    extents : list
        Dose grid extents in pixel coordintes:
        [col_pos_min, row_pos_min, col_pos_max, row_pos_max].
    dd : dict
        Dose data from dicomparser.GetDoseData.

    Returns
    -------
    list
        Dose grid extents in patient coordintes: [xmin, ymin, xmax, ymax].
    """
    return [
        dd['lut'][0][extents[0]], dd['lut'][1][extents[1]],
        dd['lut'][0][extents[2]], dd['lut'][1][extents[3]]
    ]


def get_resampled_lut(index_extents,
                      extents,
                      new_pixel_spacing,
                      min_pixel_spacing):
    """Determine the patient to pixel LUT based on new pixel spacing.

    Parameters
    ----------
    index_extents : list
        Dose grid extents as array indices.
    extents : list
        Dose grid extents in patient coordinates:
        [col_pos_min, row_pos_min, col_pos_max, row_pos_max].
    new_pixel_spacing : tuple or float
        New pixel spacing in mm (row, column).
        If float is provided, original dose grid pixel spacing must be square.
    min_pixel_spacing : tuple
        Min pixel spacing used to determine new pixel spacing (row, column).

    Returns
    -------
    tuple
        A tuple of lut lists (col lut, row lut) with the coordinates of the
        dose grid in patient coordinates

    Raises
    ------
    AttributeError
        Raised if the new pixel_spacing is not a factor of the minimum pixel
        spacing.

    Notes
    -----
    The new pixel spacing must be a factor of the original (minimum) pixel
    spacing. For example if the original pixel spacing was ``3`` mm, the new
    pixel spacing should be: ``3 / (2^n)`` mm, where ``n`` is an integer.
    This applies independently to both the row and column pixel spacing.

    If a single float value is provided it will be applied to both row and
    column. Additionally, the original dose grid pixel spacing must be square.

    Examples
    --------
    Original pixel spacing: ``3`` mm, new pixel spacing: ``0.375`` mm
    Derived via: ``(3 / 2^16) == 0.375``

    """
    if not isinstance(new_pixel_spacing, Sequence):
        if not (min_pixel_spacing[0] == min_pixel_spacing[1]):
            raise AttributeError(
                "New pixel spacing must be provided as a (row, column) tuple.")
        else:
            new_pixel_spacing = (new_pixel_spacing, new_pixel_spacing)
    if (min_pixel_spacing[0] % new_pixel_spacing[0] != 0.0):
        raise AttributeError(
            "New row pixel spacing must be a factor of %s/(2^n),"
            % min_pixel_spacing[0] +
            " where n is an integer. Value provided was %s."
            % new_pixel_spacing[0])
    if (min_pixel_spacing[1] % new_pixel_spacing[1] != 0.0):
        raise AttributeError(
            "New column pixel spacing must be a factor of %s/(2^n),"
            % min_pixel_spacing[1] +
            " where n is an integer. Value provided was %s."
            % new_pixel_spacing[1])

    # Existing number of cols, rows
    num_cols = abs(index_extents[0] - index_extents[2])
    num_rows = abs(index_extents[1] - index_extents[3])

    col_samples = round(num_cols * min_pixel_spacing[1] / new_pixel_spacing[1])
    row_samples = round(num_rows * min_pixel_spacing[0] / new_pixel_spacing[0])

    col_lut = np.linspace(
        extents[0], extents[2], int(col_samples), dtype=np.float64
    )
    row_lut = np.linspace(
        extents[1], extents[3], int(row_samples), dtype=np.float64
    )
    return col_lut, row_lut


def get_interpolated_dose(dose, z, resolution, extents):
    """Get interpolated dose for the given z, resolution & array extents.

    Parameters
    ----------
    dose : DicomParser
        A DicomParser instance of an RT Dose.
    z : float
        Index in mm of z plane of dose grid.dose
    resolution : tuple
        Interpolation resolution less than or equal to dose grid pixel spacing.
        Provided in (row, col) format.
    extents : list
        Dose grid index extents.

    Returns
    -------
    ndarray
        Interpolated dose grid with a shape larger than the input dose grid.
    """
    # Return the dose bounded by extents if interpolation is not required
    d = dose.GetDoseGrid(z)
    if not d.size:
        return d  # cannot take 2d index below if empty
    extent_dose = d[extents[1]:extents[3],
                    extents[0]:extents[2]] if len(extents) else d
    if not resolution:
        return extent_dose
    if not skimage_available:
        raise ImportError(
            "scikit-image must be installed to perform DVH interpolation.")
    scale = (np.array(dose.ds.PixelSpacing) / resolution).tolist()
    interp_dose = rescale(
        extent_dose,
        scale=scale,
        order=1,
        mode='symmetric',
        preserve_range=True,
        channel_axis=None
    )
    return interp_dose


def interpolate_between_planes(planes, n=2):
    """Interpolate n additional structure planes (segments) in between planes.

    Parameters
    ----------
    planes : dict
        RT Structure plane data from dicomparser.GetStructureCoordinates.
    n : int, optional
        Number of planes to interpolate in between the existing planes.

    Returns
    -------
    dict
        Plane data with additional keys representing interpolated planes.
    """
    keymap = {np.array([k], dtype=np.float32)[0]: k for k in planes.keys()}
    sorted_keys = np.sort(np.array(list(planes.keys()), dtype=np.float32))
    num_new_samples = (len(planes.keys()) * (n + 1)) - n
    newgrid = np.linspace(sorted_keys[0], sorted_keys[-1], num_new_samples)
    new_planes = {}
    # If the plane already exists in the dictionary, use it
    # otherwise use the closest plane
    # TODO: Add actual interpolation of structure data between planes
    for plane in newgrid:
        new_plane = sorted_keys[np.argmin(np.fabs(sorted_keys - plane))]
        new_planes[plane] = planes[keymap[new_plane]]
    return new_planes
