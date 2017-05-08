#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""unittest cases for dvhcalc."""
# test_dvhcalc.py
# Copyright (c) 2016 Aditya Panchal


import unittest
import os
from dicompylercore import dicomparser, dvhcalc
from dicompylercore.dvh import DVH
try:
    from pydicom.dataset import Dataset
    from pydicom.sequence import Sequence
except ImportError:
    from dicom.dataset import Dataset
    from dicom.sequence import Sequence
from numpy import arange

basedata_dir = "tests/testdata"
example_data = os.path.join(basedata_dir, "example_data")


class TestDVHCalc(unittest.TestCase):
    """Unit tests for DVH calculation."""

    def setUp(self):
        """Setup files for common case testing."""
        rtss_dcm = os.path.join(example_data, "rtss.dcm")
        rtdose_dcm = os.path.join(example_data, "rtdose.dcm")
        self.rtss = dicomparser.DicomParser(rtss_dcm)
        self.rtdose = dicomparser.DicomParser(rtdose_dcm)

        self.dvhs = self.rtdose.GetDVHs()

    def calc_dvh(self, key, limit=None, calculate_full_volume=True):
        """Calculate a DVH for testing."""
        # Generate the calculated DVHs
        return dvhcalc.get_dvh(
            self.rtss.ds, self.rtdose.ds, key, limit,
            calculate_full_volume=calculate_full_volume)

    def test_dvh_calculation_empty_structure_no_dose(self):
        """Test if a DVH returns an empty histogram for invalid data."""
        dvh = self.calc_dvh(2)
        self.assertEqual(dvh, DVH([0], arange(0, 2)))

    def test_dvh_calculation(self):
        """Test if cumulative DVHs can be calculated from the DICOM data."""
        dvh = self.calc_dvh(5)

        # Volume
        self.assertAlmostEqual(dvh.volume, 440.212499999)
        # Min dose bin
        self.assertAlmostEqual(dvh.bins[0], 0)
        # Max dose bin
        self.assertEqual(dvh.bins[-1], 3.100000000)
        # Max dose to structure
        self.assertAlmostEqual(dvh.max, 3.089999999)
        # Min dose to structure
        self.assertAlmostEqual(dvh.min, 0.02999999)
        # Mean dose to structure
        self.assertAlmostEqual(dvh.mean, 0.647428656)

    def test_dvh_calculation_with_dose_limit(self):
        """Test if a DVH can be calculated with a max dose limit."""
        # Set the dose limit to 500 cGy (lower than max dose)
        limitdvh = self.calc_dvh(5, limit=500)

        # Volume
        self.assertAlmostEqual(limitdvh.volume, 440.212499999)
        # Min dose bin
        self.assertAlmostEqual(limitdvh.bins[0], 0)
        # Max dose bin
        self.assertEqual(limitdvh.bins[-1], 3.100000000)
        # Max dose to structure
        self.assertAlmostEqual(limitdvh.max, 3.089999999)
        # Min dose to structure
        self.assertAlmostEqual(limitdvh.min, 0.02999999)
        # Mean dose to structure
        self.assertAlmostEqual(limitdvh.mean, 0.647428656)

        # Set the dose limit to 2000 cGy (higher than max dose)
        highlimitdvh = self.calc_dvh(5, limit=2000)
        # Max dose bin
        self.assertEqual(highlimitdvh.bins[-1], 3.100000000)

        # Set the dose limit to 1 cGy (should produce an empty histogram)
        lowlimitdvh = self.calc_dvh(5, limit=1)
        # Max dose bin
        self.assertEqual(lowlimitdvh.bins[-1], 1)

    def test_dvh_contour_outside_dose_grid(self):
        """Test if a DVH can be calculated with contours outside a dosegrid."""
        # Add a set of contours outside of the dose grid
        roi_id = 8
        roic = self.rtss.ds.ROIContourSequence[roi_id - 1]
        new_contour = Dataset()
        # Create a ContourImageSequence for the referenced Image
        new_contour.ContourImageSequence = Sequence([])
        contour_image = Dataset()
        last_contour = roic.ContourSequence[-1].ContourImageSequence[-1]
        contour_image.ReferencedSOPClassUID = \
            last_contour.ReferencedSOPClassUID
        contour_image.ReferencedSOPInstanceUID = \
            last_contour.ReferencedSOPInstanceUID
        new_contour.ContourImageSequence.append(contour_image)
        new_contour.ContourGeometricType = 'CLOSED_PLANAR'
        new_contour.NumberOfContourPoints = 4
        new_contour.ContourData = [
            0.0, -250.0, 180.0,
            5.0, -250.0, 180.0,
            5.0, -245.0, 180.0,
            0.0, -245.0, 180.0]
        roic.ContourSequence.append(new_contour)

        # Full structure volume (calculated inside/outside dose grid)
        include_vol_dvh = self.calc_dvh(8, calculate_full_volume=True)
        self.assertAlmostEqual(include_vol_dvh.volume, 0.54086538)
        # Partial volume (calculated only within dose grid)
        partial_vol_dvh = self.calc_dvh(8, calculate_full_volume=False)
        self.assertAlmostEqual(partial_vol_dvh.volume, 0.46874999)

if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
