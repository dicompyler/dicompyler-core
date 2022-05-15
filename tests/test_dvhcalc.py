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
from dicompylercore.dvhcalc import get_dvh
try:
    from pydicom.dataset import Dataset
    from pydicom.sequence import Sequence
except ImportError:
    from dicom.dataset import Dataset
    from dicom.sequence import Sequence
from numpy import arange
from numpy.testing import assert_allclose
from .util import fake_rtdose, fake_ss


basedata_dir = "tests/testdata"
example_data = os.path.join(basedata_dir, "example_data")


class TestDVHCalc(unittest.TestCase):
    """Unit tests for DVH calculation."""

    def setUp(self):
        """Set up files for common case testing."""
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

    def test_dvh_calculation_memmap(self):
        """Test if DVHs can be calculated with memmapped RT Dose."""
        dvh = dvhcalc.get_dvh(os.path.join(
            example_data, "rtss.dcm"), os.path.join(
            example_data, "rtdose.dcm"), 5, memmap_rtdose=True)
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


class TestDVHCalcDecubitus(unittest.TestCase):
    """Unit tests for DVH calculation in decubitus orientations."""

    def setUp(self):
        """Set up fake DICOM datasets used in various tests."""
        self.ss = fake_ss()
        self.dose = fake_rtdose()

    def test_nondecub(self):
        """Test that DVH is calculated correctly for standard orientation."""
        self.dose.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
        dvh = get_dvh(self.ss, self.dose, 1)
        diffl = dvh.differential
        # Counts are normalized to total, and to volume,
        # So undo that here for test dose grid.
        # 18=num dose voxels inside struct; 0.36=volume
        got_counts = diffl.counts * 18 / 0.36
        expected_counts = [0]*13 + [2, 2, 2, 0, 0, 0, 0, 0, 0, 0,
                                    2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2]
        assert_allclose(got_counts, expected_counts)

    def test_HF_decubitus_left(self):
        """Test DVH for head-first decubitus left orientation."""
        # Keep same dose grid as std orientation but pixel-spacing in X, Y same
        # For this case, use iop=[0, -1, 0, 1, 0, 0] Head first decubitus left
        # Then X = r * dr + ipp[0]
        #  and Y = -c * dc + ipp[1]
        # (https://nipy.org/nibabel/dicom/dicom_orientation.html
        # #dicom-affine-formula)
        # Change ipp y of y to new max of 19 for similar y range
        # Below show contours box of (3, 14.5) - (7, 17.5) on dose grid
        #       Y=19 18  17                  12
        # X=2   [10, 10, 10, 13, 14, 15, 16, 17],
        #               |-----------|
        #   4   [10, 10, 10, 13, 14, 15, 16, 17]
        #   6   [10, 10, 10, 13, 14, 15, 16, 17]
        #               |-----------|
        #   8   [13, 13, 13, 16, 17, 18, 19, 20]
        #  10   [14, 14, 14, 17, 18, 19, 20, 21]
        #  12   [15, 15, 15, 18, 19, 20, 21, 22]
        #  14   [16, 16, 16, 19, 20, 21, 22, 23]]

        #       Y=19 18  17                  12
        # X=2   [20, 20, 20, 23, 24, 25, 26, 27]
        #               |-----------|
        #   4   [20, 20, 20, 23, 24, 25, 26, 27]
        #   6   [20, 20, 20, 23, 24, 25, 26, 27]
        #               |-----------|
        #   8   [23, 23, 23, 26, 27, 28, 29, 30]
        #  10   [24, 24, 24, 27, 28, 29, 30, 31]
        #  12   [25, 25, 25, 28, 29, 30, 31, 32]
        #  14   [...]

        #       Y=19 18  17                  12
        # X=2   [30, 30, 30, 33, 34, 35, 36, 37]
        #               |-----------|
        #   4   [30, 30, 30, 33, 34, 35, 36, 37]
        #   6   [30, 30, 30, 33, 34, 35, 36, 37]
        #               |-----------|
        #   8   [33, 33, 33, 36, 37, 38, 39, 40]
        #   10  [34, 34, 34, 37, 38, 39, 40, 41]
        #   12  [35, 35, 35, 38, 39, 40, 41, 42]
        # X=14  [36, 36, 36, 39, 40, 41, 42, 43]

        #                          10       13 14                20
        expected_counts = [0]*10 + [2, 0, 0, 2, 2, 0, 0, 0, 0, 0, 2, 0, 0,
                                    2, 2, 0, 0, 0, 0, 0, 2, 0, 0, 2, 2]
        #                          23 24                30       33 34
        self.dose.ImagePositionPatient = [2, 19, -20]  # X Y Z top left
        self.dose.PixelSpacing = [2.0, 1.0]  # between Rows, Columns
        dvh = get_dvh(self.ss, self.dose, 1)
        diffl = dvh.differential
        # Counts are normalized to total, and to volume,
        # So undo that here for test dose grid.
        # 18=num dose voxels inside struct; 0.36=volume
        got_counts = diffl.counts * 18 / 0.36
        assert_allclose(got_counts, expected_counts)

    def test_HF_decubitus_left_structure_extents(self):
        """Test DVH for HF decubitus Lt orientation structure_extents used."""
        # Repeat test_HF_decubitus_left but with use_structure_extents
        #                          10       13 14                20
        expected_counts = [0]*10 + [2, 0, 0, 2, 2, 0, 0, 0, 0, 0, 2, 0, 0,
                                    2, 2, 0, 0, 0, 0, 0, 2, 0, 0, 2, 2]
        #                          23 24                30       33 34
        self.dose.ImagePositionPatient = [2, 19, -20]  # X Y Z top left
        self.dose.PixelSpacing = [2.0, 1.0]  # between Rows, Columns
        dvh = get_dvh(self.ss, self.dose, 1, use_structure_extents=True)
        diffl = dvh.differential
        # Counts are normalized to total, and to volume,
        # So undo that here for test dose grid.
        # 18=num dose voxels inside struct; 0.36=volume
        got_counts = diffl.counts * 18 / 0.36
        assert_allclose(got_counts, expected_counts)

    def test_HF_decubitus_right(self):
        """Test DVH for head-first decubitus right orientation."""
        # Keep same dose grid as std orientation

        self.dose.ImageOrientationPatient = [0, 1, 0, -1, 0, 0]
        self.dose.PixelSpacing = [2.0, 1.0]  # between Rows, Columns
        # original ipp = [2, 12, -20]
        # Then X = -r * dr + ipp[0], X decreases down the rows
        #  and Y = c * dc + ipp[1], Y increases across cols
        # (https://nipy.org/nibabel/dicom/dicom_orientation.html
        # #dicom-affine-formula)
        # Change ipp y of X to new max of 14 for similar y range
        self.dose.ImagePositionPatient = [14, 12, -20]  # X Y Z top left
        # Below show contours box of (3, 14.5) - (7, 17.5) on dose grid
        #       Y=12 13  14  15  16  17  18  19
        # X=14  [10, 10, 10, 13, 14, 15, 16, 17],
        #   12  [10, 10, 10, 13, 14, 15, 16, 17]
        #   10  [10, 10, 10, 13, 14, 15, 16, 17]
        #    8  [13, 13, 13, 16, 17, 18, 19, 20]
        #                   | ----------|
        #    6  [14, 14, 14, 17, 18, 19, 20, 21]
        #    4  [15, 15, 15, 18, 19, 20, 21, 22]
        #                   | ----------|
        #    2  [16, 16, 16, 19, 20, 21, 22, 23]]

        #       Y=12 13  14                  19
        # X=14  [20, 20, 20, 23, 24, 25, 26, 27]
        #   12  [20, 20, 20, 23, 24, 25, 26, 27]
        #   10  [20, 20, 20, 23, 24, 25, 26, 27]
        #    8  [23, 23, 23, 26, 27, 28, 29, 30]
        #                   | ----------|
        #    6  [24, 24, 24, 27, 28, 29, 30, 31]
        #    4  [25, 25, 25, 28, 29, 30, 31, 32]
        #                   | ----------|
        #    2  [...]

        #       Y=12 13  14                  19
        # X=14  [30, 30, 30, 33, 34, 35, 36, 37]
        #   12  [30, 30, 30, 33, 34, 35, 36, 37]
        #   10  [30, 30, 30, 33, 34, 35, 36, 37]
        #    8  [33, 33, 33, 36, 37, 38, 39, 40]
        #                   | ----------|
        #    6  [34, 34, 34, 37, 38, 39, 40, 41]
        #    4  [35, 35, 35, 38, 39, 40, 41, 42]
        #                   | ----------|
        # X= 2  [36, 36, 36, 39, 40, 41, 42, 43]

        #                           17       20
        expected_counts = [0]*17 + [1, 2, 2, 1, 0, 0, 0, 0, 0, 0,
                                    1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 1, 2, 2, 1]
        #                          27 28 29 30                   37
        dvh = get_dvh(self.ss, self.dose, 1)
        diffl = dvh.differential
        # Counts are normalized to total, and to volume,
        # So undo that here for test dose grid.
        # 18=num dose voxels inside struct; 0.36=volume
        got_counts = diffl.counts * 18 / 0.36
        assert_allclose(got_counts, expected_counts)

    def test_FF_decubitus_right(self):
        """Test DVH for feet-first decubitus right orientation."""
        self.dose.ImageOrientationPatient = [0, -1, 0, -1, 0, 0]
        self.dose.PixelSpacing = [2.0, 1.0]  # between Rows, Columns
        # original ipp = [2, 12, -20]
        # Then X = -r * dr + ipp[0], X decreases down the rows
        #  and Y = -c * dc + ipp[1], Y decreases across cols
        # (https://nipy.org/nibabel/dicom/dicom_orientation.html
        # #dicom-affine-formula)
        self.dose.ImagePositionPatient = [14, 19, 20]  # X Y Z top left
        # Below show contours box of (3, 14.5) - (7, 17.5) on dose grid
        #       Y=19 18  17  16  15  14  13  12
        # X=14  [10, 10, 10, 13, 14, 15, 16, 17],
        #   12  [10, 10, 10, 13, 14, 15, 16, 17]
        #   10  [10, 10, 10, 13, 14, 15, 16, 17]
        #    8  [13, 13, 13, 16, 17, 18, 19, 20]
        #               | ----------|
        #    6  [14, 14, 14, 17, 18, 19, 20, 21]
        #    4  [15, 15, 15, 18, 19, 20, 21, 22]
        #               | ----------|
        #    2  [16, 16, 16, 19, 20, 21, 22, 23]]

        #       Y=19 18  17  16  15  14  13  12
        # X=14  [20, 20, 20, 23, 24, 25, 26, 27]
        #   12  [20, 20, 20, 23, 24, 25, 26, 27]
        #   10  [20, 20, 20, 23, 24, 25, 26, 27]
        #    8  [23, 23, 23, 26, 27, 28, 29, 30]
        #               | ----------|
        #    6  [24, 24, 24, 27, 28, 29, 30, 31]
        #    4  [25, 25, 25, 28, 29, 30, 31, 32]
        #               | ----------|
        #    2  [...]

        #       Y=19 18  17  16  15  14  13  12
        # X=14  [30, 30, 30, 33, 34, 35, 36, 37]
        #   12  [30, 30, 30, 33, 34, 35, 36, 37]
        #   10  [30, 30, 30, 33, 34, 35, 36, 37]
        #    8  [33, 33, 33, 36, 37, 38, 39, 40]
        #               | ----------|
        #    6  [34, 34, 34, 37, 38, 39, 40, 41]
        #    4  [35, 35, 35, 38, 39, 40, 41, 42]
        #               | ----------|
        # X= 2  [36, 36, 36, 39, 40, 41, 42, 43]

        #                          14 15 16       19             24
        expected_counts = [0]*14 + [1, 1, 0, 1, 2, 1, 0, 0, 0, 0, 1, 1, 0, 1,
                                    2, 1, 0, 0, 0, 0, 1, 1, 0, 1, 2, 1]
        #                                            34
        dvh = get_dvh(self.ss, self.dose, 1)
        diffl = dvh.differential
        # Counts are normalized to total, and to volume,
        # So undo that here for test dose grid.
        # 18=num dose voxels inside struct; 0.36=volume
        got_counts = diffl.counts * 18 / 0.36
        assert_allclose(got_counts, expected_counts)

    def test_FF_decubitus_right_structure_extents(self):
        """Test DVH for FF decubitus Rt orientation using structure extents."""
        self.dose.ImageOrientationPatient = [0, -1, 0, -1, 0, 0]
        self.dose.PixelSpacing = [2.0, 1.0]  # between Rows, Columns
        self.dose.ImagePositionPatient = [14, 19, 20]  # X Y Z top left
        # see grid from test_FF_decubitus_right
        #                          14 15 16       19             24
        expected_counts = [0]*14 + [1, 1, 0, 1, 2, 1, 0, 0, 0, 0, 1, 1, 0,
                                    1, 2, 1, 0, 0, 0, 0, 1, 1, 0, 1, 2, 1]
        #                                               34
        dvh = get_dvh(self.ss, self.dose, 1, use_structure_extents=True)
        diffl = dvh.differential
        # Counts are normalized to total, and to volume,
        # So undo that here for test dose grid.
        # 18=num dose voxels inside struct; 0.36=volume
        got_counts = diffl.counts * 18 / 0.36
        assert_allclose(got_counts, expected_counts)

    def test_FF_decubitus_left(self):
        """Test DVH for feet-first decubitus left orientation."""
        self.dose.ImageOrientationPatient = [0, 1, 0, 1, 0, 0]
        self.dose.PixelSpacing = [2.0, 1.0]  # between Rows, Columns
        # original ipp = [2, 12, -20]
        # Then X = r * dr + ipp[0], X increases down the rows
        #  and Y = c * dc + ipp[1], Y increases across cols
        # (https://nipy.org/nibabel/dicom/dicom_orientation.html
        # #dicom-affine-formula)

        # In this test, we also shift Z so three structure planes use the
        #    first three dose planes rather than the middle three,
        #    just to ensure asymmetry in z direction is checked.
        #    Note, planes should really be reversed in pixel array, but doesn't
        #    matter since contour is identical on each slice.
        self.dose.ImagePositionPatient = [2, 12, 10]  # X Y Z top left
        # Below show contours box of (3, 14.5) - (7, 17.5) on dose grid
        #      Y=12  13  14  15  16  17      19
        # X=2   [ 0,  0,  0,  3,  4,  5,  6,  7],
        #                   |-----------|
        #   4   [ 0,  0,  0,  3,  4,  5,  6,  7]
        #   6   [ 0,  0,  0,  3,  4,  5,  6,  7]
        #                   |-----------|
        #   8   [ 3,  3,  3,  6,  7,  8,  9, 10]
        #  10   [ 4,  4,  4,  7,  8,  9, 10, 11]
        #  12   [ 5,  5,  5,  8,  9, 10, 11, 12]
        #  14   [ 6,  6,  6,  9, 10, 11, 12, 13]]

        #      Y=12  13  14                  19
        # X=2   [10, 10, 10, 13, 14, 15, 16, 17],
        #                   |-----------|
        #   4   [10, 10, 10, 13, 14, 15, 16, 17]
        #   6   [10, 10, 10, 13, 14, 15, 16, 17]
        #                   |-----------|
        #   8   [13, 13, 13, 16, 17, 18, 19, 20]
        #  10   [14, 14, 14, 17, 18, 19, 20, 21]
        #  12   [15, 15, 15, 18, 19, 20, 21, 22]
        #  14   [16, 16, 16, 19, 20, 21, 22, 23]]

        #      Y=12  13  14                  19
        # X=2   [20, 20, 20, 23, 24, 25, 26, 27]
        #                   |-----------|
        #   4   [20, 20, 20, 23, 24, 25, 26, 27]
        #   6   [20, 20, 20, 23, 24, 25, 26, 27]
        #                   |-----------|
        #   8   [23, 23, 23, 26, 27, 28, 29, 30]
        #  10   [24, 24, 24, 27, 28, 29, 30, 31]
        #  12   [25, 25, 25, 28, 29, 30, 31, 32]
        #  14   [...]

        #                          3
        expected_counts = [0]*3 + [2, 2, 2, 0, 0, 0, 0, 0, 0, 0,
                                   2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2]
        #                         13                            23
        dvh = get_dvh(self.ss, self.dose, 1)
        diffl = dvh.differential
        # Counts are normalized to total, and to volume,
        # So undo that here for test dose grid.
        # 18=num dose voxels inside struct; 0.36=volume
        got_counts = diffl.counts * 18 / 0.36
        assert_allclose(got_counts, expected_counts)

    def test_empty_dose_grid(self):
        """Test empty dose grid handled correctly."""
        # See #274, prior to fixes this raised IndexError from
        #  get_interpolated_dose() getting empty array from GetDoseGrid()
        # Use z value to force no dose grid at that value
        #  Otherwise make like decub example
        self.dose.ImagePositionPatient = [2, 19, -1020]  # X Y Z top left
        self.dose.PixelSpacing = [2.0, 1.0]  # between Rows, Columns

        # 1 = roi number
        dvh = get_dvh(self.ss, self.dose, 1, use_structure_extents=True)
        self.assertTrue('Empty DVH' in dvh.notes)

    def test_not_implemented_orientations(self):
        """Test unhandled orientations raise NotImplementedError."""
        self.dose.ImageOrientationPatient = [0.7071, 0.7071, 0, 1, 0, 0]
        with self.assertRaises(NotImplementedError):
            _ = get_dvh(self.ss, self.dose, 1)


if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
