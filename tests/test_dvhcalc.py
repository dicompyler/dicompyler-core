#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""unittest cases for dvhcalc."""
# test_dvhcalc.py
# Copyright (c) 2016 Aditya Panchal


import unittest
import os
from dicompylercore import dicomparser, dvhcalc
from dicompylercore.dvh import DVH
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

    def calc_dvh(self, key, limit=None):
        """Calculate a DVH for testing."""
        # Generate the calculated DVHs
        return dvhcalc.get_dvh(self.rtss.ds, self.rtdose.ds, key, limit)

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
        # Set the dose limit to 100 cGy
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

if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
