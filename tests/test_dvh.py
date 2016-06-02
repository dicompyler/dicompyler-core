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
        rtdose = dicomparser.DicomParser(filename=rtdose_dcm)

        cls.dvh = dvh.DVH.from_dicom_dvh(rtdose.ds, 7)
        cls.rx_dose = 14

    def test_raw_data_dvh(self):
        """Test if a DVH can be created from raw data with a [0, 1] bin."""
        self.assertEqual(dvh.DVH.from_data(1, 1), dvh.DVH([1], [1]))
        self.assertEqual(
            dvh.DVH.from_data(1, 1).__repr__(),
            "DVH(cumulative, 1 bins: [0:1] Gy, volume: 0 cm3)")

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

    def test_dvh_properties(self):
        """Test if the DVH properties can be derived."""
        self.assertEqual(self.dvh.max, 14.569999999999734)
        self.assertEqual(self.dvh.min, 14.069999999999744)
        self.assertEqual(self.dvh.mean, 14.285830178442305)
        self.assertEqual(self.dvh.volume, 12.809180549338702)

    @unittest.skipUnless(mpl_available, "Matplotlib not installed")
    def test_plotting(self):
        """Test if the DVH can be plotted."""
        self.assertEqual(self.dvh.plot(), self.dvh)

if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
