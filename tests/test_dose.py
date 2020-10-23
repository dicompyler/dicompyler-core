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
    from pydicom import read_file as read_dicom
except ImportError:
    from dicom.dataset import Dataset
    from dicom.sequence import Sequence
    from dicom import read_file as read_dicom
from numpy import arange, zeros
from numpy.testing import (
    assert_array_almost_equal,
    assert_array_equal,
    assert_raises,
)
import warnings
from dicompylercore.config import mpl_available, scipy_available

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
        self.rtdose_dcm = os.path.join(example_data, "rtdose.dcm")
        self.rtdose = dicomparser.DicomParser(self.rtdose_dcm)
        self.dosegrid = dose.DoseGrid(self.rtdose_dcm)

    def test_modality_check(self):
        """Test non-RTDOSE raises AttributeError"""
        ds = dicomparser.DicomParser(self.rtdose_dcm).ds
        ds.Modality = "RTPLAN"
        with self.assertRaises(AttributeError):
            dose.DoseGrid(ds)

    def test_shape(self):
        """Test dose grid shape extraction from the DICOM data."""
        self.assertEqual(self.dosegrid.shape, (194, 129, 98))

    def test_scale(self):
        """Test dose grid resolution extraction from the DICOM data."""
        assert_array_equal(self.dosegrid.scale, [2.5, 2.5, 3.0])

    def test_offset(self):
        """Test dose grid resolution extraction from the DICOM data."""
        assert_array_equal(
            self.dosegrid.offset, [-228.6541915, -419.2444776, -122.4407]
        )

    def test_is_coincident(self):
        """Test spatial coincidence of two dose grids"""

        # Self coincidence
        other = dose.DoseGrid(self.rtdose_dcm)
        self.assertTrue(self.dosegrid.is_coincident(other))

        # ImagePositionPatient (offset)
        other = dose.DoseGrid(self.rtdose_dcm)
        other.ds.ImagePositionPatient[0] += 1
        self.assertFalse(self.dosegrid.is_coincident(other))

        # PixelSpacing
        other = dose.DoseGrid(self.rtdose_dcm)
        other.ds.PixelSpacing[0] += 1
        self.assertFalse(self.dosegrid.is_coincident(other))

        # GridFrameOffsetVector
        check = arange(0, 98) * 3
        assert_array_equal(self.dosegrid.ds.GridFrameOffsetVector, check)
        other = dose.DoseGrid(self.rtdose_dcm)
        other.ds.GridFrameOffsetVector[0] += 1
        self.assertFalse(self.dosegrid.is_coincident(other))

    def test_add_overload(self):
        """Test the overloaded __add__ operator"""
        other = dose.DoseGrid(self.rtdose_dcm)
        dose_sum = self.dosegrid + other
        assert_array_equal(self.dosegrid.dose_grid, other.dose_grid)
        assert_array_equal(dose_sum.dose_grid, self.dosegrid.dose_grid * 2)

    def test_multiply_overload(self):
        """Test the overloaded __mul__ operator"""
        scaled_dosegrid = self.dosegrid * 2
        assert_array_equal(
            scaled_dosegrid.dose_grid, self.dosegrid.dose_grid * 2
        )

    def test_right_multiply_overload(self):
        """Test the overloaded __rmul__ operator"""
        scaled_dosegrid = 2 * self.dosegrid
        assert_array_equal(
            scaled_dosegrid.dose_grid, self.dosegrid.dose_grid * 2
        )

    def test_multiply(self):
        """Directly test the multiply function"""
        dosegrid = dose.DoseGrid(self.rtdose_dcm)
        dosegrid.multiply(2)
        assert_array_equal(dosegrid.dose_grid, self.dosegrid.dose_grid * 2)

        # Check that a negative factor raises NotImplementedError
        with self.assertRaises(NotImplementedError):
            dosegrid.multiply(-1)

    def test_dose_direct_sum(self):
        """Test the direct summation of two coincident dose grids"""
        other = dose.DoseGrid(self.rtdose_dcm)
        other._direct_sum(self.dosegrid)
        assert_array_equal(other.dose_grid, self.dosegrid.dose_grid * 2)

    @unittest.skipUnless(scipy_available, "scipy not installed")
    def test_dose_interp_sum(self):
        """Test the interpolated summation of two coincident dose grids"""
        other = dose.DoseGrid(self.rtdose_dcm)
        other._interp_sum(self.dosegrid)
        assert_array_almost_equal(other.dose_grid, self.dosegrid.dose_grid * 2)

    @unittest.skipUnless(scipy_available, "scipy not installed")
    def test_interp_entire_grid(self):
        """Test interp_entire_grid of two non-coincident dose grids"""
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

    @unittest.skipUnless(scipy_available, "scipy not installed")
    def test_interp_param(self):
        """Test interp summation of two non-coincident dose grids with
        non-default interp parameters"""
        # Interp Sum equality, entire grid in one operation
        other_ds = read_dicom(self.rtdose_dcm)
        other_ds.ImagePositionPatient[0] += 0.0000005
        other = dose.DoseGrid(other_ds)
        dosegrid = dose.DoseGrid(self.rtdose_dcm, order=2, mode="nearest")
        dosegrid.add(other)
        assert_array_not_equal(dosegrid.dose_grid, self.dosegrid.dose_grid * 2)
        assert_array_almost_equal(
            dosegrid.dose_grid, self.dosegrid.dose_grid * 2
        )

        # Interp Sum inequality, entire grid in one operation
        other_ds.ImagePositionPatient[0] += 0.0000005
        other = dose.DoseGrid(other_ds)
        dosegrid = dose.DoseGrid(self.rtdose_dcm, order=2, mode="nearest")
        dosegrid.add(other)
        assert_array_not_equal(dosegrid.dose_grid, self.dosegrid.dose_grid * 2)
        assert_array_not_almost_equal(
            dosegrid.dose_grid, self.dosegrid.dose_grid * 2
        )

        self.assertEqual(dosegrid.interp_param["order"], 2)
        self.assertEqual(dosegrid.interp_param["mode"], "nearest")

    def test_add_attr_mismatch(self):
        """Test add fails with mismatched DoseSummationType"""
        warnings.filterwarnings("ignore")
        other = dose.DoseGrid(self.rtdose_dcm)
        other.ds.DoseSummationType = "%s1" % other.ds.DoseSummationType

        with self.assertRaises(NotImplementedError):
            other.add(self.dosegrid)
        warnings.filterwarnings("default")

    @unittest.skipUnless(mpl_available, "Matplotlib not installed")
    def test_show(self):
        """Test if the dose grid can be shown."""
        import matplotlib.pyplot as plt

        plt.ion()
        self.assertEqual(self.dosegrid.show().ds, self.dosegrid.ds)

    def test_set_dicom_tag_value(self):
        """Test set_dicom_tag_value by tag and keyword"""
        # Edit existing tag by keyword
        ds = dose.DoseGrid(self.rtdose_dcm).ds
        self.assertNotEqual(str(ds.PatientID), "DoseTestByKeyword")
        dose.set_dicom_tag_value(ds, "PatientID", "DoseTestByKeyword")
        self.assertEqual(str(ds.PatientID), "DoseTestByKeyword")

        # Edit existing tag by tag
        self.assertNotEqual(str(ds.PatientID), "DoseTestByTag")
        dose.set_dicom_tag_value(ds, 0x00100020, "DoseTestByTag")
        self.assertEqual(str(ds.PatientID), "DoseTestByTag")

        # Create a new Tag
        self.assertFalse(hasattr(ds, "PatientComments"))
        dose.set_dicom_tag_value(ds, "PatientComments", "CommentsTest")
        self.assertEqual(str(ds.PatientComments), "CommentsTest")

    def test_add_dicom_sequence(self):
        """Test add_dicom_sequence by appending or creating a new sequence"""

        # Add new sequence
        ds = dose.DoseGrid(self.rtdose_dcm).ds
        self.assertFalse(hasattr(ds, "ReferencedInstanceSequence"))
        seq_data = {"ReferencedSOPClassUID": "TestUID1"}
        dose.add_dicom_sequence(ds, "ReferencedInstanceSequence", seq_data)
        self.assertEqual(
            str(ds.ReferencedInstanceSequence[0].ReferencedSOPClassUID),
            "TestUID1",
        )

        # Append to existing sequence
        seq_data = {"ReferencedSOPClassUID": "TestUID2"}
        dose.add_dicom_sequence(ds, "ReferencedInstanceSequence", seq_data)
        self.assertEqual(
            str(ds.ReferencedInstanceSequence[1].ReferencedSOPClassUID),
            "TestUID2",
        )

    def test_validate_attr_equality(self):
        """Test validate_attr_equality"""

        # Check equivalence of test attr of two TestObj objects
        obj_1 = TestObj("test value")
        obj_2 = TestObj("test value")
        self.assertTrue(dose.validate_attr_equality(obj_1, obj_2, "test"))

        # Check test attr of two TestObj objects are not equal
        obj_2.test = "test fail"
        warnings.filterwarnings("ignore")
        self.assertFalse(dose.validate_attr_equality(obj_1, obj_2, "test"))
        warnings.filterwarnings("default")

    def test_save_dcm(self):
        """Test save DoseGrid to DICOM"""

        dosegrid = dose.DoseGrid(self.rtdose_dcm)
        self.assertFalse(hasattr(dosegrid.ds, "ContentDate"))
        self.assertFalse(hasattr(dosegrid.ds, "ContentTime"))

        dosegrid2 = dose.DoseGrid(self.rtdose_dcm)
        dosegrid.add(dosegrid2)  # ensure other_sop_class_uid is set

        filepath = os.path.join(example_data, "dose_write_test.dcm")
        dosegrid.save_dcm(filepath)

        dosegrid_new = dose.DoseGrid(filepath)  # load new dosegrid from file
        self.assertTrue(hasattr(dosegrid_new.ds, "ContentDate"))
        self.assertTrue(hasattr(dosegrid_new.ds, "ContentTime"))

        os.remove(filepath)

    def test_boundary_dose(self):
        """Check boundary dose calculations"""
        self.assertEqual(self.dosegrid.max_boundary_dose, 0.138544)
        self.assertAlmostEqual(
            self.dosegrid.max_boundary_relative_dose, 0.009437111038635319
        )

    def test_non_uniform_dose_grid_scale(self):
        """Check that a non-uniform dose grid is detected"""
        ds = dose.DoseGrid(self.rtdose_dcm).ds
        ds.GridFrameOffsetVector[0] += 1
        dosegrid = dose.DoseGrid(ds)

        with self.assertRaises(NotImplementedError):
            dosegrid.scale

    def test_max_boundary_value(self):
        """Test the max_boundary_value function"""
        arr = zeros([3, 3, 3])
        arr[1, 1, 1] = 10

        arr[0, 1, 1] = 1
        self.assertEqual(dose.max_boundary_value(arr), 1)
        arr[-1, 1, 1] = 2
        self.assertEqual(dose.max_boundary_value(arr), 2)

        arr[1, 1, 0] = 3
        self.assertEqual(dose.max_boundary_value(arr), 3)
        arr[1, 1, -1] = 4
        self.assertEqual(dose.max_boundary_value(arr), 4)

        arr[1, 0, 1] = 5
        self.assertEqual(dose.max_boundary_value(arr), 5)
        arr[1, -1, 1] = 6
        self.assertEqual(dose.max_boundary_value(arr), 6)


class TestObj(object):
    """Simple test class for test_validate_attr_equality"""

    def __init__(self, init_value):
        self.test = init_value


if __name__ == "__main__":
    import sys

    sys.exit(unittest.main())
