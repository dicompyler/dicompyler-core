#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""unittest cases for dvhcalc."""
# test_dvhcalc.py
# Copyright (c) 2016 Aditya Panchal


import unittest
import os
from dicompylercore import dicomparser, dvhcalc
from numpy import array
from numpy.testing import assert_array_equal, assert_array_almost_equal

basedata_dir = "tests/testdata"
example_data = os.path.join(basedata_dir, "example_data")

# import logging
# import logging.handlers
# logger = logging.getLogger('dicompylercore.dvhcalc')
# logger.setLevel(logging.DEBUG)
# ch = logging.StreamHandler()
# ch.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
# logger.addHandler(ch)


class TestDVHCalc(unittest.TestCase):
    """Unit tests for DVH calculation."""

    def setUp(self):
        """Setup files for common case testing."""
        rtss_dcm = os.path.join(example_data, "rtss.dcm")
        rtdose_dcm = os.path.join(example_data, "rtdose.dcm")
        self.rtss = dicomparser.DicomParser(filename=rtss_dcm)
        self.rtdose = dicomparser.DicomParser(filename=rtdose_dcm)

        self.structures = self.rtss.GetStructures()
        self.dvhs = self.rtdose.GetDVHs()

    def test_dvh_calculation(self):
        """Test if DVHs can be calculated from the DICOM data."""

        # Generate the calculated DVHs
        key = 5
        structure = self.structures[key]
        structure['planes'] = self.rtss.GetStructureCoordinates(key)
        structure['thickness'] = self.rtss.CalculatePlaneThickness(
            structure['planes'])
        dvh = dvhcalc.get_dvh(structure, self.rtdose)

        self.assertAlmostEqual(dvh['data'][0], 440.2312499999)
        self.assertAlmostEqual(dvh['data'][-1], 0.018749999999)


if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
