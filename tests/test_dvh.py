#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""unittest cases for dvh."""
# test_dvh.py
# Copyright (c) 2016 Aditya Panchal


import unittest
import os
from dicompylercore import dvh, dicomparser
from numpy import array, arange
from numpy.testing import assert_array_equal

mpl_available = True
try:
    import matplotlib.pyplot as plt
except:
    mpl_available = False
else:
    plt.ioff()

basedata_dir = "tests/testdata"
example_data = os.path.join(basedata_dir, "example_data")


class TestDVH(unittest.TestCase):
    """Unit tests for the DVH module."""

    @classmethod
    def setUpClass(cls):
        """Setup files for common case testing."""
        rtdose_dcm = os.path.join(example_data, "rtdose.dcm")
        rtdose = dicomparser.DicomParser(rtdose_dcm)

        cls.dvh = dvh.DVH.from_dicom_dvh(rtdose.ds, 9)
        cls.rx_dose = 14

    def test_raw_data_dvh(self):
        """Test if a DVH can be created from raw data."""
        self.assertEqual(dvh.DVH.from_data(1, 1), dvh.DVH([1], [1]))
        self.assertEqual(
            repr(dvh.DVH.from_data(1, 1)),
            "DVH(cumulative, 1 bins: [0:1] Gy, volume: 1 cm3, "
            "name: None, rx_dose: 0 Gy)")
        assert_array_equal(dvh.DVH.from_data(0, 1).bins, array([0, 0]))
        assert_array_equal(dvh.DVH.from_data(5, 2).bins, array([0, 2, 4, 5]))

    def test_raw_data_dvh_max_bins(self):
        """Test if a DVH can be created from raw data with [0, 5] bin."""
        max_dvh = dvh.DVH.from_data([0, 5])
        assert_array_equal(max_dvh.counts, array([1, 0, 0, 0, 1]))
        assert_array_equal(max_dvh.bins, arange(0, 6))

    def test_differential_dvh(self):
        """Test if a cumulative DVH can be converted to a differential DVH."""
        self.assertAlmostEqual(
            self.dvh.counts.max(), self.dvh.differential.counts.sum())
        self.assertEqual(
            self.dvh.differential, self.dvh.differential.differential)

    def test_cumulative_dvh(self):
        """Test if a differential DVH can be converted to a cumulative DVH."""
        self.assertEqual(
            self.dvh.cumulative, self.dvh.differential.cumulative)

    def test_absolute_relative_dose_dvh(self):
        """Test if an absolute and relative dose DVH can be generated."""
        self.assertEqual(
            self.dvh.absolute_dose(self.rx_dose),
            self.dvh.relative_dose(self.rx_dose).absolute_dose(self.rx_dose))
        self.assertEqual(
            self.dvh.relative_dose(self.rx_dose).relative_dose(self.rx_dose),
            self.dvh.absolute_dose(self.rx_dose).relative_dose(self.rx_dose))

    def test_absolute_relative_volume_dvh(self):
        """Test if an absolute and relative volume DVH can be generated."""
        self.assertEqual(
            self.dvh.absolute_volume(self.dvh.volume),
            self.dvh.relative_volume.absolute_volume(self.dvh.volume))
        self.assertEqual(
            self.dvh.relative_volume.relative_volume,
            self.dvh.absolute_volume(self.dvh.volume).relative_volume)

    def test_absolute_relative_full_conversion(self):
        """Test if an abs / relative volume / dose DVH can be generated."""
        a = self.dvh
        b = a.differential.relative_volume.absolute_volume(a.volume).cumulative
        c = a.relative_volume.differential.absolute_volume(a.volume).cumulative
        d = a.relative_volume.absolute_volume(a.volume).differential.cumulative
        e = a.differential.relative_volume.cumulative.absolute_volume(a.volume)
        self.assertEqual(a, b)
        self.assertEqual(b, c)
        self.assertEqual(c, d)
        self.assertEqual(d, e)

        rx = self.rx_dose
        f = b.relative_dose(rx).absolute_dose(rx).differential.relative_volume
        g = b.relative_dose(rx).relative_volume.differential.absolute_dose(rx)
        h = b.relative_volume.differential.relative_dose(rx).absolute_dose(rx)
        i = b.differential.relative_dose(rx).relative_volume.absolute_dose(rx)
        self.assertEqual(f, g)
        self.assertEqual(g, h)
        self.assertEqual(h, i)

        # Test if rx_dose is included in initial constructor
        i.rx_dose = rx
        j = i.relative_dose().absolute_dose().differential.relative_volume
        k = i.relative_dose().relative_volume.differential.absolute_dose()
        l = i.relative_volume.differential.relative_dose().absolute_dose()
        m = i.differential.relative_dose().relative_volume.absolute_dose()
        self.assertEqual(i, j)
        self.assertEqual(j, k)
        self.assertEqual(k, l)
        self.assertEqual(l, m)

        # Test if rx_dose is not included in initial constructor
        # but is accessed as if it is provided
        with self.assertRaises(AttributeError):
            self.dvh.relative_dose()
        with self.assertRaises(AttributeError):
            self.dvh.relative_dose(14).absolute_dose()

    def test_dvh_dosescaling(self):
        """Test if the DVH Dose Scaling is applied properly."""
        rtdose_dcm = os.path.join(example_data, "rtdose.dcm")
        rtdose = dicomparser.DicomParser(rtdose_dcm)
        rtdose.ds.DVHSequence[7].DVHDoseScaling = 2
        scaleddvh = dvh.DVH.from_dicom_dvh(rtdose.ds, 9)
        self.assertEqual(self.dvh.max * 2, scaleddvh.max)
        self.assertEqual(self.dvh.min * 2, scaleddvh.min)
        self.assertEqual(self.dvh.mean * 2, scaleddvh.mean)
        self.assertEqual(self.dvh.volume, scaleddvh.volume)

    def test_dvh_properties(self):
        """Test if the DVH properties can be derived."""
        self.assertEqual(self.dvh.max, 14.579999999999734)
        self.assertEqual(self.dvh.min, 14.069999999999744)
        self.assertEqual(self.dvh.mean, 14.285830178442307)
        self.assertEqual(self.dvh.volume, 12.809180549338803)

    def test_dvh_value(self):
        """Test if the DVHValue class works as expected."""
        self.assertEqual(str(dvh.DVHValue(100)), '100.00')
        self.assertEqual(str(dvh.DVHValue(100, 'Gy')), '100.00 Gy')
        self.assertEqual(
            repr(dvh.DVHValue(100, 'Gy')),
            "dvh.DVHValue(100, 'Gy')")

    def test_dvh_statistics(self):
        """Test if the DVH statistics can be calculated."""
        self.dvh.rx_dose = self.rx_dose
        self.assertEqual(
            self.dvh.volume_constraint(0),
            dvh.DVHValue(12.809180549338601, 'cm3'))
        self.assertEqual(
            self.dvh.volume_constraint(100),
            dvh.DVHValue(12.809180549338601, 'cm3'))
        self.assertEqual(
            self.dvh.volume_constraint(105),
            dvh.DVHValue(0.0, 'cm3'))
        self.assertEqual(
            self.dvh.volume_constraint(14, 'Gy'),
            dvh.DVHValue(12.809180549338601, 'cm3'))
        self.assertEqual(
            self.dvh.volume_constraint(100, 'Gy'),
            dvh.DVHValue(0.0, 'cm3'))
        self.assertEqual(
            self.dvh.dose_constraint(100),
            dvh.DVHValue(14.059999999999745, 'Gy'))
        self.assertEqual(
            self.dvh.dose_constraint(90),
            dvh.DVHValue(14.169999999999742, 'Gy'))
        self.assertEqual(
            self.dvh.dose_constraint(0.02, 'cc'),
            dvh.DVHValue(14.529999999999735, 'Gy'))
        self.assertEqual(
            self.dvh.dose_constraint(15, 'cc'),
            dvh.DVHValue(0.0, 'Gy'))

    def test_dvh_statistics_shorthand(self):
        """Test if the DVH statistics can be accessed via shorthand."""
        self.assertEqual(
            self.dvh.v100, dvh.DVHValue(12.809180549338601, 'cm3'))
        self.assertEqual(
            self.dvh.v14Gy, dvh.DVHValue(12.809180549338601, 'cm3'))
        self.assertEqual(
            self.dvh.D90, dvh.DVHValue(14.169999999999742, 'Gy'))
        self.assertEqual(
            self.dvh.d2cc, dvh.DVHValue(14.389999999999738, 'Gy'))

    def test_dvh_statistics_shorthand_fail(self):
        """Test if the DVH statistics shorthand fail on invalid syntaxes."""
        with self.assertRaises(AttributeError):
            self.dvh.v100agy

    def test_dvh_describe(self):
        """Test if the DVH statistics summary can be generated."""
        self.assertEqual(self.dvh.describe(), None)
        self.assertEqual(self.dvh.relative_dose(self.rx_dose).describe(), None)

    def test_dvh_compare(self):
        """Test if the DVH comparsion summary can be generated."""
        self.dvh.name = "test"
        self.assertEqual(self.dvh.compare(self.dvh), None)
        self.assertEqual(self.dvh.relative_dose(
            self.rx_dose).compare(
                self.dvh.relative_dose(self.rx_dose)), None)
        with self.assertRaises(AttributeError):
            self.dvh.relative_dose(self.rx_dose).compare(self.dvh)

    @unittest.skipUnless(mpl_available, "Matplotlib not installed")
    def test_plotting(self):
        """Test if the DVH can be plotted."""
        self.assertEqual(self.dvh.plot(), self.dvh)
        self.dvh.name = "test"
        self.assertEqual(self.dvh.plot(), self.dvh)

    def test_dvh_statistics_with_no_counts(self):
        subject = dvh.DVH(array([]), array([0]))
        self.assertEqual(subject.max, 0)
        self.assertEqual(subject.min, 0)
        self.assertEqual(subject.mean, 0)

    def test_dose_constraint_with_no_counts(self):
        subject = dvh.DVH(array([]), array([0]))
        subject.dose_constraint(1)

    def test_dvh_statistics_with_zero_volume(self):
        subject = dvh.DVH(array([0, 0]), array([0, 1]))
        self.assertEqual(subject.max, 0)
        self.assertEqual(subject.min, 0)
        self.assertEqual(subject.mean, 0)

    def test_dose_constraint_with_zero_volume(self):
        subject = dvh.DVH(array([0, 0]), array([0, 1]))
        subject.dose_constraint(1)


if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
