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
#    https://github.com/cutright/DVH-Analytics/

from copy import deepcopy
import numpy as np
from dicompylercore import dicomparser
from pydicom.uid import generate_uid
from pydicom.datadict import dictionary_VR, keyword_dict
from dicompylercore.config import (
    dicompyler_uid_prefix_rtdose,
    mpl_available,
    scipy_available,
)
from datetime import datetime
from pydicom.sequence import Sequence
from pydicom.dataset import Dataset
from warnings import warn

if scipy_available:
    from scipy.ndimage import map_coordinates


class DoseGrid:
    """Class that stores DICOM-RT dose grids, performs addition/scaling."""

    def __init__(
        self,
        rt_dose,
        order=1,
        mode="constant",
        cval=0.0,
    ):
        """ Initialization of a DoseGrid from a DICOM-RT Dose file or dataset.

        Parameters
        ----------
        rt_dose : pydicom Dataset or filename
            DICOM RT Dose used to determine the structure dose grid data.
        order : int, optional
            The order of the spline interpolation (if needed), default is 1.
            The order has to be in the range 0-5.
            0: the nearest grid point, 1: trilinear, 2 to 5: spline
            See scipy.ndimage.map_coordinates documentation for more details
        mode : 'constant' or 'nearest', optional
            The mode parameter determines how the other dose grid is extended
            beyond its boundaries. Default is ``'constant'``. Behavior for
            these values is as follows:

            ``'constant'`` (k k k k | a b c d | k k k k)
                The other dose grid is extended by filling all values beyond
                the edge with the same constant value, defined by the cval
                parameter.
            ``'nearest'`` (a a a a | a b c d | d d d d)
                The input is extended by replicating the last pixel.

            Additional modes are available, see scipy.ndimage.map_coordinates
            documentation for more details.
        cval : scalar, optional
            Value to fill past edges of input if mode is ‘constant’.
            Default is 0.0.
        """

        self.ds = dicomparser.DicomParser(rt_dose).ds

        self.interp_param = {"order": order, "mode": mode, "cval": cval}

        self.summation_type = None
        self.sop_class_uid = getattr(self.ds, 'SOPClassUID', '')
        self.sop_instance_uid = getattr(self.ds, 'SOPInstanceUID', '')
        self.other_sop_class_uid = None
        self.other_sop_instance_uid = None

        if self.ds.Modality == "RTDOSE":
            self.x_axis = (
                np.arange(self.ds.Columns) * self.ds.PixelSpacing[0]
                + self.ds.ImagePositionPatient[0]
            )
            self.y_axis = (
                np.arange(self.ds.Rows) * self.ds.PixelSpacing[1]
                + self.ds.ImagePositionPatient[1]
            )
            self.z_axis = (
                np.array(self.ds.GridFrameOffsetVector)
                + self.ds.ImagePositionPatient[2]
            )

            # x and z are swapped in the pixel_array
            pixel_array = self.ds.pixel_array * self.ds.DoseGridScaling
            self.dose_grid = np.swapaxes(pixel_array, 0, 2)
        else:
            raise AttributeError(
                "The DoseGrid class requires an RTDOSE file or dataset. "
                "%s was detected" % self.ds.Modality
            )

    ####################################################
    # Basic properties
    ####################################################
    @property
    def shape(self):
        """Get the x, y, z dimensions of the dose grid"""
        return (
            self.ds.Columns,
            self.ds.Rows,
            len(self.ds.GridFrameOffsetVector),
        )

    @property
    def axes(self):
        """Get the x, y, z axes of the dose grid (in mm)"""
        return [self.x_axis, self.y_axis, self.z_axis]

    @property
    def scale(self):
        """Get the dose grid resolution (xyz)"""
        diffs = np.diff(self.ds.GridFrameOffsetVector)
        if not np.all(np.isclose(diffs, [diffs[0]]*len(diffs))):
            raise NotImplementedError(
                "Non-uniform GridFrameOffsetVector detected. Interpolated "
                "summation of non-uniform dose-grid scales is not supported."
            )
        return np.array(
            [
                self.ds.PixelSpacing[0],
                self.ds.PixelSpacing[1],
                self.ds.GridFrameOffsetVector[1]
                - self.ds.GridFrameOffsetVector[0],
            ]
        )

    @property
    def offset(self):
        """Get the coordinates of the dose grid origin (mm)"""
        return np.array(self.ds.ImagePositionPatient, dtype="float")

    @property
    def max_boundary_dose(self):
        """Get the max boundary dose"""
        return max_boundary_value(self.dose_grid)

    @property
    def max_boundary_relative_dose(self):
        return self.max_boundary_dose / np.max(self.dose_grid)

    ####################################################
    # Tools
    ####################################################
    def __add__(self, other):
        """Overload + operator to sum this dose grid with the other dose grid

        Parameters
        ----------
        other : DoseGrid
            Another DoseGrid object.
        """
        new = deepcopy(self)
        new.add(other)
        return new

    def __mul__(self, factor):
        """Overload * operator to scale this dose grid by the provided factor

        Parameters
        ----------
        factor : int, float
            Scale the dose grid by this value.
        """
        new = deepcopy(self)
        new.multiply(factor)
        return new

    def __rmul__(self, factor):
        return self.__mul__(factor)

    def multiply(self, factor):
        """
        Scale the dose grid.

        Parameters
        ----------
        factor : int, float
            Multiply the dose grid by this factor.
        """

        if factor < 0:
            raise NotImplementedError("Negative doses are not supported.")

        self.dose_grid *= factor
        self.dose_grid_post_processing()

    def dose_grid_post_processing(self, other=None):
        """Set the pixel data and store UIDs from other DoseGrid"""
        self.set_pixel_data()
        if hasattr(self.ds, "DVHSequence"):
            del self.ds.DVHSequence
        if other is not None:
            self.other_sop_class_uid = other.sop_class_uid
            self.other_sop_instance_uid = other.sop_instance_uid

    def is_coincident(self, other):
        """Check dose grid spatial coincidence.

        Parameters
        ----------
        other : DoseGrid
            Another DoseGrid object.
        """
        return (
            self.ds.PixelSpacing == other.ds.PixelSpacing
            and self.ds.ImagePositionPatient == other.ds.ImagePositionPatient
            and self.ds.pixel_array.shape == other.ds.pixel_array.shape
            and self.ds.GridFrameOffsetVector == other.ds.GridFrameOffsetVector
        )

    def set_pixel_data(self):
        """Update the PixelData with the current dose_grid"""
        self.ds.BitsAllocated = 32
        self.ds.BitsStored = 32
        self.ds.HighBit = 31
        self.ds.DoseGridScaling = (
            np.max(self.dose_grid) / np.iinfo(np.uint32).max
        )
        pixel_data = (
            np.swapaxes(self.dose_grid, 0, 2) / self.ds.DoseGridScaling
        )
        self.ds.PixelData = np.uint32(pixel_data).tobytes()

    def save_dcm(self, file_path):
        """Save the pydicom.FileDataset to file"""
        self.update_dicom_tags()
        self.ds.save_as(file_path)

    def get_ijk_points(self, other_axes):
        """Convert axes from another DoseGrid into ijk of this DoseGrid.

        Parameters
        ----------
        other_axes : list
            The x, y, and z axis arrays.

        Returns
        -------
        np.vstack
            Array of other_axes in this ijk space.
        """
        ijk_axes = [
            (np.array(axis) - self.offset[a]) / self.scale[a]
            for a, axis in enumerate(other_axes)
        ]
        j, i, k = np.meshgrid(ijk_axes[1], ijk_axes[0], ijk_axes[2])
        return np.vstack((i.ravel(), j.ravel(), k.ravel()))

    ####################################################
    # Dose Summation
    ####################################################
    def add(self, other, force=False):
        """
        Add another dose grid to this dose grid, with interpolation if needed

        Parameters
        ----------
        other : DoseGrid
            Another DoseGrid object.
        force : bool
            Set to True to ignore differences in DoseSummationType, DoseType,
            DoseUnits, ImageOrientationPatient
        """

        attrs = [
            "DoseSummationType",
            "DoseType",
            "DoseUnits",
            "ImageOrientationPatient",
        ]
        attr_check = [
            validate_attr_equality(self.ds, other.ds, attr) for attr in attrs
        ]
        if not force and not all(attr_check):
            mismatches = [
                attr for i, attr in enumerate(attrs) if attr_check[i]
            ]
            raise NotImplementedError(
                "Dose summation of dose grids with these mismatched "
                "attributes is not recommended: %s. Use "
                "DoseGrid.add(other, force=True) to ignore"
                % ",".join(mismatches)
            )

        if self.is_coincident(other):
            self._direct_sum(other)
        else:
            if not scipy_available:
                raise ImportError(
                    "scipy must be installed to perform interpolated dose sum."
                )
            self._interp_sum(other)

    def _direct_sum(self, other):
        """Directly sum two coincident dose grids

        Parameters
        ----------
        other: DoseGrid
            Another DoseGrid object.
        """
        self.dose_grid += other.dose_grid
        self.summation_type = "DIRECT"

        self.dose_grid_post_processing(other)

    def _interp_sum(self, other):
        """
        Interpolate the other dose grid to this dose grid's axes,
        then perform direct summation

        Parameters
        ----------
        other: DoseGrid
            Another DoseGrid object.
        """

        self.dose_grid += self.interp_entire_grid(other)
        self.summation_type = "INTERPOLATED"

        self.dose_grid_post_processing(other)

    def interp_entire_grid(self, other):
        """
        Interpolate the other dose grid to this dose grid's axes in one
        operation

        Parameters
        ----------
        other: DoseGrid
            Another DoseGrid object.

        Returns
        -------
        np.array
            The other dose grid interpolated to this dose grid's axes
        """
        return map_coordinates(
            input=other.dose_grid,
            coordinates=other.get_ijk_points(self.axes),
            **self.interp_param
        ).reshape(self.shape)

    def update_dicom_tags(self):
        """Update DICOM UIDs, Content Date/Time, and Dose Comment"""

        # Store the source SOPClassUID and SOPInstanceUID
        seq_data = {
            "ReferencedSOPClassUID": self.sop_class_uid,
            "ReferencedSOPInstanceUID": self.sop_instance_uid,
        }
        add_dicom_sequence(self.ds, "ReferencedInstanceSequence", seq_data)

        if self.other_sop_class_uid is not None:
            seq_data = {
                "ReferencedSOPClassUID": self.other_sop_class_uid,
                "ReferencedSOPInstanceUID": self.other_sop_instance_uid,
            }
            add_dicom_sequence(self.ds, "ReferencedInstanceSequence", seq_data)

        # Create a new SOPInstanceUID
        set_dicom_tag_value(
            self.ds,
            "SOPInstanceUID",
            generate_uid(prefix=dicompyler_uid_prefix_rtdose),
        )

        # Store the dose summation type in the DoseComment tag
        if self.summation_type:
            set_dicom_tag_value(
                self.ds, "DoseComment", "%s SUMMATION" % self.summation_type
            )

        # Update the Date and Time tags
        now = datetime.now()
        set_dicom_tag_value(self.ds, "ContentDate", now.strftime("%Y%m%d"))
        set_dicom_tag_value(self.ds, "ContentTime", now.strftime("%H%M%S"))

    def show(self, z=None):
        """Show the dose grid using Matplotlib if present.

        Parameters
        ----------
        z : float, optional
            slice position to display initially, by default None

        """
        if not mpl_available:
            raise ImportError(
                "Matplotlib could not be loaded. Install and try again.")
            return self
        import matplotlib.pyplot as plt
        from matplotlib.widgets import Slider

        # Extract the list of planes (z) from the dose grid
        planes = (
            np.array(self.ds.GridFrameOffsetVector)
            * self.ds.ImageOrientationPatient[0]
            * self.ds.ImageOrientationPatient[4]
        ) + self.ds.ImagePositionPatient[2]

        # Set up the plot
        fig = plt.figure()
        ax = fig.add_subplot(111)
        rtdose = dicomparser.DicomParser(self.ds)

        # Get the middle slice if the z is not provided
        z = planes[planes.size // 2] if z is None else z
        zplane = rtdose.GetDoseGrid(z) * self.ds.DoseGridScaling
        # Flag to invert slider min/max if GFOV is decreasing (i.e. FFS)
        reverse = planes[0] > planes[-1]
        im = ax.imshow(zplane, cmap="jet",)

        # Create a slider to change the (z)
        axslice = fig.add_axes([0.34, 0.01, 0.50, 0.02])
        slider = Slider(
            ax=axslice,
            label="Slice Position (mm):",
            valmin=planes[-1] if reverse else planes[0],
            valmax=planes[0] if reverse else planes[-1],
            valinit=z,
            valstep=np.diff(planes)[0],
        )

        def updateslice(z):
            """Update the data to show on the plot."""
            im.set_data(rtdose.GetDoseGrid(z) * self.ds.DoseGridScaling)
            plt.draw()

        slider.on_changed(updateslice)
        plt.show()
        return self


def set_dicom_tag_value(ds, tag, value):
    """Set or update a DICOM tag value in the pydicom dataset.

    Parameters
    ----------
    ds : pydicom Dataset
        The pydicom dataset for the tag to be added/updated to.
    tag : str, int or tuple
        DICOM tag or keyword to be added.
    value : any
        New value for the tag's element.
    """
    try:
        ds[tag].value = value
    except KeyError:
        if tag in keyword_dict:  # Keyword provided rather than int or tuple
            tag = keyword_dict[tag]
        ds.add_new(tag, dictionary_VR(tag), value)


def add_dicom_sequence(ds, seq_keyword, data_set_dict):
    """Add a sequence to a data set.

    Parameters
    ----------
    ds : pydicom Dataset
        The pydicom dataset for the sequence to be added to.
    seq_keyword : str
        The DICOM keyword for the sequence.
    data_set_dict : dict
        Dictionary of tags and values for the sequence element.
    """
    seq_ds = Dataset()
    for tag, value in data_set_dict.items():
        set_dicom_tag_value(seq_ds, tag, value)

    if hasattr(ds, seq_keyword):
        getattr(ds, seq_keyword).append(seq_ds)
    else:
        setattr(ds, seq_keyword, Sequence([seq_ds]))


def validate_attr_equality(obj_1, obj_2, attr):
    """Assess the equality of the provided attr between two objects.
    Send warning if unequal.

    Parameters
    ----------
    obj_1 : object
        Any object with an `attr` attribute that is comparable by !=
    obj_2 : object
        Any object with an `attr` attribute that is comparable by !=
    attr : str
        The attribute to be compared between obj_1 and obj_2
    """
    val_1 = getattr(obj_1, attr)
    val_2 = getattr(obj_2, attr)
    if val_1 != val_2:
        warn("Different %s values detected:\n%s\n%s" % (attr, val_1, val_2))
        return False
    return True


def max_boundary_value(arr):
    """Get the greatest value on the boundary of a 3D numpy array

    Parameters
    ----------
    arr : numpy.array
        Any 3-dimensional array-like object

    Returns
    -------
    float
        Maximum value along any side of the input array
    """
    return np.max(
        [
            np.max([np.max(arr[i, :, :]) for i in [0, -1]]),
            np.max([np.max(arr[:, j, :]) for j in [0, -1]]),
            np.max([np.max(arr[:, :, k]) for k in [0, -1]]),
        ]
    )
