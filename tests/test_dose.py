#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""unittest cases for dose."""
# test_dvhdose.py
# Copyright (c) 2016-2020 Aditya Panchal
# Copyright (c) 2020 Dan Cutright


from __future__ import division
import unittest
import os
from dicompylercore import dicomparser, dvhcalc, dose
from dicompylercore.dvh import DVH
from dicompylercore.config import scipy_available
try:
    from pydicom.dataset import Dataset
    from pydicom.sequence import Sequence
except ImportError:
    from dicom.dataset import Dataset
    from dicom.sequence import Sequence
from numpy import arange, zeros
from numpy.testing import assert_array_almost_equal, assert_array_equal

basedata_dir = "tests/testdata"
example_data = os.path.join(basedata_dir, "example_data")


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
        """Test if the dose grid shape can be extracted from the DICOM data."""
        shape = self.dosegrid.shape

        # x-dimension
        self.assertEqual(shape[0], 194)
        # y-dimension
        self.assertEqual(shape[1], 129)
        # z-dimension
        self.assertEqual(shape[2], 98)

    def test_scale(self):
        """Test if the dose grid resolution can be extracted from the DICOM data."""
        scale = self.dosegrid.scale

        # x-dimension
        self.assertEqual(scale[0], 2.5)
        # y-dimension
        self.assertEqual(scale[1], 2.5)
        # z-dimension
        self.assertEqual(scale[2], 3.)

    def test_offset(self):
        """Test if the dose grid resolution can be extracted from the DICOM data."""
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
        other = dose.DoseGrid(self.rtdose_dcm)
        other.ds.GridFrameOffsetVector[0] += 1
        self.assertEqual(self.dosegrid.is_coincident(other), False)

        # GridFrameOffsetVector
        # other = dose.DoseGrid(self.rtdose_dcm)
        # other.ds.pixel_array = zeros([10, 10, 10])  # AttributeError: can't set attribute
        # self.assertEqual(self.dosegrid.is_coincident(other), False)

    def test_dose_direct_sum(self):
        # Direct Sum with a factor
        other = dose.DoseGrid(self.rtdose_dcm)
        other.direct_sum(self.dosegrid, other_factor=2)
        assert_array_equal(other.dose_grid, self.dosegrid.dose_grid * 3)

    def test_dose_interp_sum(self):
        # Interp Sum with a factor
        other = dose.DoseGrid(self.rtdose_dcm)
        other.interp_sum(self.dosegrid, other_factor=2)
        assert_array_almost_equal(other.dose_grid, self.dosegrid.dose_grid * 3)

    def test_interp_entire_grid(self):
        # Interp Sum, entire grid in one operation
        other = dose.DoseGrid(self.rtdose_dcm)
        other.interp_entire_grid(self.dosegrid)
        assert_array_almost_equal(other.dose_grid, self.dosegrid.dose_grid)

    def test_interp_by_block(self):
        # Interp Sum by block
        other = dose.DoseGrid(self.rtdose_dcm)
        other.interp_by_block(self.dosegrid)
        assert_array_almost_equal(other.dose_grid, self.dosegrid.dose_grid)


if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
