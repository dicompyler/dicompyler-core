#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""unittest cases for dose."""
# test_dvhdose.py
# Copyright (c) 2016-2020 Aditya Panchal
# Copyright (c) 2020 Dan Cutright


from __future__ import division
import unittest
import os
from dicompylercore import dicomparser, dose
try:
    from pydicom.dataset import Dataset
    from pydicom.sequence import Sequence
except ImportError:
    from dicom.dataset import Dataset
    from dicom.sequence import Sequence
from numpy import arange
from numpy.testing import (
    assert_array_almost_equal,
    assert_array_equal,
    assert_raises,
)

basedata_dir = "tests/testdata"
example_data = os.path.join(basedata_dir, "example_data")


def assert_array_not_equal(arr_1, arr_2):
    return assert_raises(AssertionError, assert_array_equal, arr_1, arr_2)


def assert_array_not_almost_equal(arr_1, arr_2):
    return assert_raises(
        AssertionError, assert_array_almost_equal, arr_1, arr_2
    )


class TestDose(unittest.TestCase):
    """Unit tests for Dose calculations."""

    def setUp(self):
        """Setup files for common case testing."""
        rtss_dcm = os.path.join(example_data, "rtss.dcm")
        self.rtdose_dcm = os.path.join(example_data, "rtdose.dcm")
        self.rtss = dicomparser.DicomParser(rtss_dcm)
        self.rtdose = dicomparser.DicomParser(self.rtdose_dcm)

        self.dvhs = self.rtdose.GetDVHs()

        self.dosegrid = dose.DoseGrid(self.rtdose_dcm)

    def test_shape(self):
        """Test dose grid shape extraction from the DICOM data."""
        shape = self.dosegrid.shape

        # x-dimension
        self.assertEqual(shape[0], 194)
        # y-dimension
        self.assertEqual(shape[1], 129)
        # z-dimension
        self.assertEqual(shape[2], 98)

    def test_scale(self):
        """Test dose grid resolution extraction from the DICOM data."""
        scale = self.dosegrid.scale

        # x-dimension
        self.assertEqual(scale[0], 2.5)
        # y-dimension
        self.assertEqual(scale[1], 2.5)
        # z-dimension
        self.assertEqual(scale[2], 3.0)

    def test_offset(self):
        """Test dose grid resolution extraction from the DICOM data."""
        offset = self.dosegrid.offset

        # x-dimension
        self.assertEqual(offset[0], -228.6541915)
        # y-dimension
        self.assertEqual(offset[1], -419.2444776)
        # z-dimension
        self.assertEqual(offset[2], -122.4407)

    def test_is_coincident(self):
        # Self coincidence
        other = dose.DoseGrid(self.rtdose_dcm)
        self.assertEqual(self.dosegrid.is_coincident(other), True)

        # ImagePositionPatient (offset)
        other = dose.DoseGrid(self.rtdose_dcm)
        other.ds.ImagePositionPatient[0] += 1
        self.assertEqual(self.dosegrid.is_coincident(other), False)

        # PixelSpacing
        other = dose.DoseGrid(self.rtdose_dcm)
        other.ds.PixelSpacing[0] += 1
        self.assertEqual(self.dosegrid.is_coincident(other), False)

        # GridFrameOffsetVector
        check = arange(0, 98) * 3
        assert_array_equal(self.dosegrid.ds.GridFrameOffsetVector, check)
        other = dose.DoseGrid(self.rtdose_dcm)
        other.ds.GridFrameOffsetVector[0] += 1
        self.assertEqual(self.dosegrid.is_coincident(other), False)

    def test_add_overload(self):
        # Overloaded __add__ operator
        other = dose.DoseGrid(self.rtdose_dcm)
        dose_sum = self.dosegrid + other
        assert_array_equal(self.dosegrid.dose_grid, other.dose_grid)
        assert_array_equal(dose_sum.dose_grid, self.dosegrid.dose_grid * 2)

    def test_multiply_overload(self):
        # Overloaded __mul__ operator
        scaled_dosegrid = self.dosegrid * 2
        assert_array_equal(
            scaled_dosegrid.dose_grid, self.dosegrid.dose_grid * 2
        )

    def test_right_multiply_overload(self):
        # Overloaded __rmul__ operator
        scaled_dosegrid = 2 * self.dosegrid
        assert_array_equal(
            scaled_dosegrid.dose_grid, self.dosegrid.dose_grid * 2
        )

    def test_dose_direct_sum(self):
        # Direct Sum with a factor
        other = dose.DoseGrid(self.rtdose_dcm)
        other.direct_sum(self.dosegrid)
        assert_array_equal(other.dose_grid, self.dosegrid.dose_grid * 2)

    def test_dose_interp_sum(self):
        # Interp Sum with a factor
        other = dose.DoseGrid(self.rtdose_dcm)
        other.interp_sum(self.dosegrid)
        assert_array_almost_equal(other.dose_grid, self.dosegrid.dose_grid * 2)

    def test_interp_entire_grid(self):
        # Interp Sum equality, entire grid in one operation
        other = dose.DoseGrid(self.rtdose_dcm)
        other.x_axis += 0.0000005  # perturb to ensure interpolation is used
        interp_grid_1 = other.interp_entire_grid(self.dosegrid)
        assert_array_not_equal(interp_grid_1, self.dosegrid.dose_grid)
        assert_array_almost_equal(interp_grid_1, self.dosegrid.dose_grid)

        # Interp Sum inequality, entire grid in one operation
        other.x_axis += 0.0000005  # should cause assert_array_equal failure
        interp_grid_2 = other.interp_entire_grid(self.dosegrid)
        assert_array_not_equal(interp_grid_2, self.dosegrid.dose_grid)
        assert_array_not_almost_equal(interp_grid_2, self.dosegrid.dose_grid)

    def test_add_attr_mismatch(self):
        other = dose.DoseGrid(self.rtdose_dcm)
        other.ds.DoseSummationType = ""
        failure_detected = False
        try:
            other.add(self.dosegrid)
        except NotImplementedError:
            failure_detected = True
        self.assertEqual(failure_detected, True)


if __name__ == "__main__":
    import sys

    sys.exit(unittest.main())
