#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""unittest cases for dvhcalc."""
# test_dvhcalc.py
# Copyright (c) 2016-2018 Aditya Panchal


from __future__ import division
import unittest
import os
from dicompylercore import dicomparser, dvhcalc
from dicompylercore.config import skimage_available
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

    def calc_dvh(self, key, limit=None,
                 calculate_full_volume=True,
                 use_structure_extents=False,
                 interpolation_resolution=None,
                 interpolation_segments=0):
        """Calculate a DVH for testing."""
        # Generate the calculated DVHs
        dvh = dvhcalc.get_dvh(
            self.rtss.ds, self.rtdose.ds, key, limit,
            calculate_full_volume=calculate_full_volume,
            use_structure_extents=use_structure_extents,
            interpolation_resolution=interpolation_resolution,
            interpolation_segments_between_planes=interpolation_segments)
        dvh.dose_units = 'Gy'
        return dvh

    def create_new_contour(self, roi_id, extents, z):
        """Create a new contour sequence for the given ROI id."""

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
        xmin, ymin, xmax, ymax = extents
        new_contour.ContourData = [
            xmin, ymin, z,
            xmax, ymin, z,
            xmax, ymax, z,
            xmin, ymax, z
        ]
        roic.ContourSequence.append(new_contour)

    def test_dvh_calculation_empty_structure_no_dose(self):
        """Test if a DVH returns an empty histogram for invalid data."""
        dvh = self.calc_dvh(2)
        self.assertEqual(dvh, DVH([0], arange(0, 2)))

    def test_dvh_calculation(self):
        """Test if cumulative DVHs can be calculated from the DICOM data."""
        dvh = self.calc_dvh(5)

        # Volume
        self.assertAlmostEqual(dvh.volume, 440.23124999)
        # Min dose bin
        self.assertAlmostEqual(dvh.bins[0], 0)
        # Max dose bin
        self.assertEqual(dvh.bins[-1], 3.1)
        # Max dose to structure
        self.assertAlmostEqual(dvh.max, 3.1)
        # Min dose to structure
        self.assertAlmostEqual(dvh.min, 0.03)
        # Mean dose to structure
        self.assertAlmostEqual(dvh.mean, 0.6475329)

    def test_dvh_calculation_with_dose_limit(self):
        """Test if a DVH can be calculated with a max dose limit."""
        # Set the dose limit to 500 cGy (lower than max dose)
        limitdvh = self.calc_dvh(5, limit=500)

        # Volume
        self.assertAlmostEqual(limitdvh.volume, 440.23124999)
        # Min dose bin
        self.assertAlmostEqual(limitdvh.bins[0], 0)
        # Max dose bin
        self.assertEqual(limitdvh.bins[-1], 3.1)
        # Max dose to structure
        self.assertAlmostEqual(limitdvh.max, 3.1)
        # Min dose to structure
        self.assertAlmostEqual(limitdvh.min, 0.03)
        # Mean dose to structure
        self.assertAlmostEqual(limitdvh.mean, 0.6475329)

        # Set the dose limit to 2000 cGy (higher than max dose)
        highlimitdvh = self.calc_dvh(5, limit=2000)
        # Max dose bin
        self.assertEqual(highlimitdvh.bins[-1], 3.1)

        # Set the dose limit to 1 cGy (should produce an empty histogram)
        lowlimitdvh = self.calc_dvh(5, limit=1)
        # Max dose bin
        self.assertEqual(lowlimitdvh.bins[-1], 1)

    def test_dvh_contour_outside_dose_grid(self):
        """Test if a DVH can be calculated with contours outside a dosegrid."""
        # Add a set of contours outside of the dose grid
        self.create_new_contour(8, [0.0, -250.0, 5.0, -245.0], 180.0)

        # Full structure volume (calculated inside/outside dose grid)
        include_vol_dvh = self.calc_dvh(8, calculate_full_volume=True)
        self.assertAlmostEqual(include_vol_dvh.volume, 0.56249999)
        # Partial volume (calculated only within dose grid)
        partial_vol_dvh = self.calc_dvh(8, calculate_full_volume=False)
        self.assertAlmostEqual(partial_vol_dvh.volume, 0.48749999)

    @unittest.skipUnless(skimage_available, "scikit-image not installed")
    def test_dvh_with_in_plane_interpolation(self):
        """Test if DVH can be calculated using in plane interpolation."""
        interp_dvh = self.calc_dvh(
            8, use_structure_extents=True,
            interpolation_resolution=(2.5 / 8))

        # Volume
        self.assertAlmostEqual(interp_dvh.volume, 0.51590551)
        # Min dose bin
        self.assertAlmostEqual(interp_dvh.bins[0], 0)
        # Max dose bin
        self.assertEqual(interp_dvh.bins[-1], 12.98)
        # Max dose to structure
        self.assertAlmostEqual(interp_dvh.max, 12.98)
        # Min dose to structure
        self.assertAlmostEqual(interp_dvh.min, 1.32)
        # Mean dose to structure
        self.assertAlmostEqual(interp_dvh.mean, 7.695116550116536)

    @unittest.skipUnless(skimage_available, "scikit-image not installed")
    def test_dvh_with_in_plane_interpolation_non_square_pixel_spacing(self):
        """Test non-square pixel spacing DVH calculation with interpolation."""
        interp_dvh = self.calc_dvh(
            8, use_structure_extents=True,
            interpolation_resolution=((2.5 / 8), (2.5 / 16)))

        # Volume
        self.assertAlmostEqual(interp_dvh.volume, 0.51215152)
        # Min dose bin
        self.assertAlmostEqual(interp_dvh.bins[0], 0)
        # Max dose bin
        self.assertEqual(interp_dvh.bins[-1], 13.01)
        # Max dose to structure
        self.assertAlmostEqual(interp_dvh.max, 13.01)
        # Min dose to structure
        self.assertAlmostEqual(interp_dvh.min, 1.37)
        # Mean dose to structure
        self.assertAlmostEqual(interp_dvh.mean, 7.660532286212908)

        # Fake irregular pixel spacing to test resampled LUT errors
        # for non square pixel spacing
        print(self.rtdose.ds.PixelSpacing)
        self.rtdose.ds.PixelSpacing = [2.0, 3.0]

        # Test that a non-sequence resolution is invalid
        # for non-square pixel spacing
        with self.assertRaises(AttributeError):
            self.calc_dvh(
                8, use_structure_extents=True,
                interpolation_resolution=(2.5 / 8))

        # Test row incorrect new pixel spacing
        with self.assertRaises(AttributeError):
            self.calc_dvh(
                8, use_structure_extents=True,
                interpolation_resolution=((2.1 / 8), (3.0 / 16)))

        # Test column incorrect pixel spacing
        with self.assertRaises(AttributeError):
            self.calc_dvh(
                8, use_structure_extents=True,
                interpolation_resolution=((2.0 / 8), (3.1 / 8)))

    def test_dvh_with_structure_extents(self):
        """Test if DVH calculation is same as normal with structure extents."""
        orig_dvh = self.calc_dvh(8)
        structure_extents_dvh = self.calc_dvh(8, use_structure_extents=True)
        self.assertEqual(orig_dvh, structure_extents_dvh)

    def test_dvh_with_structure_extents_larger_than_dose_grid(self):
        """Test DVH calculation using large structure structure extents."""
        # Add a set of contours larger than the dose grid plane
        self.create_new_contour(3, [-230.0, -520.0, 260.0, 0.0], 24.56)

        structure_extents_dvh = self.calc_dvh(3, use_structure_extents=True)
        self.assertAlmostEqual(structure_extents_dvh.volume, 464.40000)

    def test_dvh_with_in_plane_interpolation_sampling_fail(self):
        """Test if DVH calculation fails when the sampling rate is invalid."""
        with self.assertRaises(AttributeError):
            self.calc_dvh(
                8, use_structure_extents=False,
                interpolation_resolution=(3 / 8))

    def test_dvh_calculation_with_interpolation_between_planes(self):
        """Test if DVH can be calculated using interpolation between planes."""
        dvh = self.calc_dvh(8, interpolation_segments=2)

        # Volume
        self.assertAlmostEqual(dvh.volume, 0.47499999)
        # Min dose bin
        self.assertAlmostEqual(dvh.bins[0], 0)
        # Max dose bin
        self.assertEqual(dvh.bins[-1], 10.0)
        # Max dose to structure
        self.assertAlmostEqual(dvh.max, 10.0)
        # Min dose to structure
        self.assertAlmostEqual(dvh.min, 2.03)
        # Mean dose to structure
        self.assertAlmostEqual(dvh.mean, 6.4767105)


if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
