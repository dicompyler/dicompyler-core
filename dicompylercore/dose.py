#!/usr/bin/env python
# -*- coding: utf-8 -*-
# dose.py
"""Routines to access and modify DICOM RT Dose."""
# Copyright (c) 2009-2016 Aditya Panchal
# Copyright (c) 2019-2020 Dan Cutright
# This file is part of dicompyler-core, released under a BSD license.
#    See the file license.txt included with this distribution, also
#    available at https://github.com/dicompyler/dicompyler-core/
#
# This code was adapted from dicom_dose_sum.py from DVH Analytics:
#    https://github.com/cutright/DVH-Analytics/blob/master/dvha/tools/dicom_dose_sum.py

from copy import deepcopy
import numpy as np
from scipy.ndimage import map_coordinates
from dicompylercore import dicomparser
from pydicom.uid import generate_uid
from pydicom.datadict import dictionary_VR, keyword_dict
from dicompylercore.config import dicompyler_uid_prefix_rtdose
from datetime import datetime
from pydicom.sequence import Sequence
from pydicom.dataset import Dataset


class DoseGrid:
    """
    Class to easily access commonly used attributes of a DICOM dose grid and perform make modifications

    Example: Add two dose grids
        grid_1 = DoseGrid(dose_file_1)
        grid_2 = DoseGrid(dose_file_2)
        grid_sum = grid_1 + grid_2
        grid_sum.save_dcm(some_file_path)

    """
    def __init__(self, rt_dose, order=1, try_full_interp=True, interp_block_size=50000):
        """
        :param rt_dose: pydicom Dose Dataset or filename to a DICOM RT-Dose
        :param order: The order of the spline interpolation. The order has to be in the range 0-5.
                      0: the nearest grid point, 1: trilinear, 2 to 5: spline
        :param try_full_interp: If true, will attempt to interpolate the entire grid at once before calculating one
                                block at a time (block size defined in self.interp_by_block)
        :type try_full_interp: bool
        :param interp_block_size: calculate this many points at a time if not try_full_interp or MemoryError
        :type interp_block_size: int
        """

        self.ds = dicomparser.DicomParser(rt_dose).ds
        self.order = order
        self.try_full_interp = try_full_interp
        self.interp_block_size = interp_block_size

        self.summation_type = None
        self.sop_class_uid = self.ds.SOPClassUID
        self.sop_instance_uid = self.ds.SOPInstanceUID
        self.other_sop_class_uid = None
        self.other_sop_instance_uid = None

        if self.ds:
            self.__set_axes()

    def __set_axes(self):
        self.x_axis = np.arange(self.ds.Columns) * self.ds.PixelSpacing[0] + self.ds.ImagePositionPatient[0]
        self.y_axis = np.arange(self.ds.Rows) * self.ds.PixelSpacing[1] + self.ds.ImagePositionPatient[1]
        self.z_axis = np.array(self.ds.GridFrameOffsetVector) + self.ds.ImagePositionPatient[2]

        # x and z are swapped in the pixel_array
        self.dose_grid = np.swapaxes(self.ds.pixel_array * self.ds.DoseGridScaling, 0, 2)

    ####################################################
    # Basic properties
    ####################################################
    @property
    def shape(self):
        """Get the x, y, z dimensions of the dose grid"""
        return tuple([self.ds.Columns, self.ds.Rows, len(self.ds.GridFrameOffsetVector)])

    @property
    def axes(self):
        """Get the x, y, z axes of the dose grid (in mm)"""
        return [self.x_axis, self.y_axis, self.z_axis]

    @property
    def scale(self):
        """Get the dose grid resolution (xyz)"""
        return np.array([self.ds.PixelSpacing[0],
                         self.ds.PixelSpacing[1],
                         self.ds.GridFrameOffsetVector[1] - self.ds.GridFrameOffsetVector[0]])

    @property
    def offset(self):
        """Get the coordinates of the dose grid origin (mm)"""
        return np.array(self.ds.ImagePositionPatient, dtype='float')

    @property
    def points(self):
        """Get all of the points in the dose grid"""
        y, x, z = np.meshgrid(self.y_axis, self.x_axis, self.z_axis)
        points = np.vstack((x.ravel(), y.ravel(), z.ravel()))
        return points.transpose()

    ####################################################
    # Tools
    ####################################################
    def __add__(self, other):
        """Add other dose grid"""
        new = deepcopy(self)
        new.add(other)
        return new

    def __mul__(self, factor):
        """Multiply this dose grid by a factor"""
        new = deepcopy(self)
        new.scale_by_factor(factor)
        return new

    def __rmul__(self, factor):
        return self.__mul__(factor)

    def is_coincident(self, other):
        """Check dose grid coincidence, if True a direct summation is appropriate"""
        return self.ds.ImagePositionPatient == other.ds.ImagePositionPatient and \
               self.ds.pixel_array.shape == other.ds.pixel_array.shape and \
               self.ds.PixelSpacing == other.ds.PixelSpacing and \
               self.ds.GridFrameOffsetVector == other.ds.GridFrameOffsetVector

    def set_pixel_data(self):
        """Update the PixelData in the pydicom.FileDataset with the current self.dose_grid"""
        self.ds.BitsAllocated = 32
        self.ds.BitsStored = 32
        self.ds.HighBit = 31
        self.ds.DoseGridScaling = np.max(self.dose_grid) / np.iinfo(np.uint32).max
        pixel_data = np.swapaxes(self.dose_grid, 0, 2) / self.ds.DoseGridScaling
        self.ds.PixelData = np.uint32(pixel_data).tostring()

    def save_dcm(self, file_path):
        """Save the pydicom.FileDataset to file"""
        self.ds.save_as(file_path)

    def get_ijk_points(self, other_axes):
        """
        Convert axes from another DoseGrid into ijk of this DoseGrid
        :param other_axes: the x, y, and z axis arrays
        :type other_axes: list
        :return: np.vstack of other_axes in this ijk space
        """
        ijk_axes = [(np.array(axis) - self.offset[a]) / self.scale[a] for a, axis in enumerate(other_axes)]
        j, i, k = np.meshgrid(ijk_axes[1], ijk_axes[0], ijk_axes[2])
        return np.vstack((i.ravel(), j.ravel(), k.ravel()))

    def scale_by_factor(self, factor):
        """
        Scale the dose grid
        :param factor: scale the dose grid by this factor
        :type factor: int or float
        """
        self.dose_grid *= factor
        self.summation_post_processing()

    ####################################################
    # Dose Summation
    ####################################################
    def add(self, other, other_factor=1):
        """
        Add another 3D dose grid to this 3D dose grid, with interpolation if needed
        :param other: another DoseGrid
        :type other: DoseGrid
        :param other_factor
        """
        if self.is_coincident(other):
            self.direct_sum(other, other_factor)
        else:
            self.interp_sum(other, other_factor)

    def direct_sum(self, other, other_factor=1):
        """Directly sum two dose grids (only works if both are coincident)"""
        self.dose_grid += other.dose_grid * other_factor
        self.summation_type = "DIRECT"

        self.summation_post_processing(other)

    def interp_sum(self, other, other_factor=1):
        """
        Interpolate the other dose grid to this dose grid's axes, then directly sum
        :param other: another DoseGrid
        :type other: DoseGrid
        :param other_factor:
        """
        other_grid = None
        if self.try_full_interp:
            try:
                other_grid = self.interp_entire_grid(other)
            except MemoryError:
                pass
        if other_grid is None:
            other_grid = self.interp_by_block(other)

        self.dose_grid += other_grid * other_factor
        self.summation_type = "INTERPOLATED"

        self.summation_post_processing(other)

    def summation_post_processing(self, other=None):
        self.set_pixel_data()
        if other is not None:
            self.other_sop_class_uid = other.sop_class_uid
            self.other_sop_instance_uid = other.sop_instance_uid
        self.update_dicom_tags()

    def interp_entire_grid(self, other):
        """
        Interpolate the other dose grid to this dose grid's axes
        :param other: another DoseGrid
        :type other: DoseGrid
        """
        points = other.get_ijk_points(self.axes)
        return map_coordinates(input=other.dose_grid, coordinates=points, order=self.order).reshape(self.shape)

    def interp_by_block(self, other):
        """
        Interpolate the other dose grid to this dose grid's axes, calculating one block at a time
        The block is defined at the init of this class, default is 50,000 points at a time
        :param other: another DoseGrid
        :type other: DoseGrid
        """
        points = other.get_ijk_points(self.axes)
        point_count = np.product(self.shape)
        other_grid = np.zeros(point_count)

        block_count = int(np.floor(point_count / self.interp_block_size))

        for i in range(block_count):
            start = i * self.interp_block_size
            end = (i+1) * self.interp_block_size if i + 1 < block_count else -1
            other_grid[start:end] = map_coordinates(input=other.dose_grid, coordinates=points[start:end],
                                                    order=self.order)

        return other_grid.reshape(self.shape)

    def update_dicom_tags(self):

        # Store the source SOPClassUID and SOPInstanceUID
        seq_data = {'ReferencedSOPClassUID': self.sop_class_uid,
                    'ReferencedSOPInstanceUID': self.sop_instance_uid}
        add_dicom_sequence(self.ds, 'ReferencedInstanceSequence', seq_data)

        if self.other_sop_class_uid is not None:
            seq_data = {'ReferencedSOPClassUID': self.other_sop_class_uid,
                        'ReferencedSOPInstanceUID': self.other_sop_instance_uid}
            add_dicom_sequence(self.ds, 'ReferencedInstanceSequence', seq_data)

        # Create a new SOPInstanceUID
        set_dicom_tag_value(self.ds, 'SOPInstanceUID', generate_uid(prefix=dicompyler_uid_prefix_rtdose))

        # Store the dose summation type in the DoseComment tag
        if self.summation_type:
            set_dicom_tag_value(self.ds, 'DoseComment', "%s SUMMATION" % self.summation_type)

        # Update the Date and Time tags
        now = datetime.now()
        set_dicom_tag_value(self.ds, 'ContentDate', now.strftime("%Y%m%d"))
        set_dicom_tag_value(self.ds, 'ContentTime', now.strftime("%H%M%S"))


def set_dicom_tag_value(ds, tag, value):
    """
    Edit a DICOM tag value, create a new element if it does not exist in the pydicom dataset
    :param ds: pydicom dataset
    :param tag: DICOM tag or keyword
    :type tag: int, tuple, or str
    :param value: new value for the tag's element
    """
    try:
        ds[tag].value = value
    except KeyError:
        if tag in keyword_dict:  # Keyword provided rather than int or tuple
            tag = keyword_dict[tag]
        try:
            ds.add_new(tag, dictionary_VR(tag), value)
        except Exception as e:
            print("Invalid DICOM tag/keyword: %s\n%s" % (tag, e))


def add_dicom_sequence(ds, seq_keyword, data_set_dict):
    """
    Add a sequence to a data set
    :param ds: pydicom dataset
    :param seq_keyword: the DICOM keyword for the sequence
    :param data_set_dict: dictionary of tags and values for the sequence element
    """

    seq_ds = Dataset()
    for tag, value in data_set_dict.items():
        set_dicom_tag_value(seq_ds, tag, value)

    if hasattr(ds, seq_keyword):
        getattr(ds, seq_keyword).append(seq_ds)
    else:
        setattr(ds, seq_keyword, Sequence([seq_ds]))



