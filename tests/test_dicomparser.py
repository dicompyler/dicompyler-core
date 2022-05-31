#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""unittest cases for dicomparser."""
# test_database.py
# Copyright (c) 2016 Aditya Panchal


import unittest
import os
from dicompylercore import dicomparser
from dicompylercore.config import pil_available, shapely_available
try:
    from pydicom.multival import MultiValue as mv
    from pydicom.valuerep import DSfloat
except ImportError:
    from dicom.multival import MultiValue as mv
    from dicom.valuerep import DSfloat
from numpy import array, arange
from numpy.testing import assert_array_equal, assert_array_almost_equal

basedata_dir = "tests/testdata"
example_data = os.path.join(basedata_dir, "example_data")
nonimage_data = os.path.join(basedata_dir, "ecg")
charsettests = os.path.join(basedata_dir, "charsettests")

# import logging
# import logging.handlers
# logger = logging.getLogger('dicompylercore.dicomparser')
# logger.setLevel(logging.DEBUG)
# ch = logging.StreamHandler()
# ch.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
# logger.addHandler(ch)


class TestCommon(unittest.TestCase):
    """Unit tests for common cases."""

    def setUp(self):
        """Setup files for common case testing."""
        ct_0_dcm = os.path.join(example_data, "ct.0.dcm")
        self.dp = dicomparser.DicomParser(ct_0_dcm)

    def test_file_import(self):
        """Test if a standard DICOM file can be parsed."""
        self.assertEqual(self.dp.GetSOPClassUID(), 'ct')

    def test_missing_file_meta(self):
        """Test if a standard DICOM file can be parsed."""
        file_meta_original = self.dp.ds.file_meta.copy()
        self.dp.ds.file_meta.clear()
        dp = dicomparser.DicomParser(self.dp.ds)
        file_meta_fixed = dp.ds.file_meta
        self.assertEqual(file_meta_original, file_meta_fixed)

    def test_dataset_import(self):
        """Test if a pydicom dataset file can be parsed."""
        dp1 = dicomparser.DicomParser(self.dp.ds)
        self.assertEqual(self.dp.ds, dp1.ds)

    def test_study_info(self):
        """Test if the study info can be parsed."""
        study = {
            'description': 'No description',
            'date': '19010101',
            'id': '2.16.840.1.113662.2.12.0.3057.1241703565.35',
            'time': '000000'
        }
        self.assertEqual(self.dp.GetStudyInfo(), study)

    def test_series_info(self):
        """Test if the series info can be parsed."""
        series = {
            'date': '19010101',
            'description': 'No description',
            'id': '2.16.840.1.113662.2.12.0.3057.1241703565.43',
            'study': '2.16.840.1.113662.2.12.0.3057.1241703565.35',
            'referenceframe': '2.16.840.1.113662.2.12.0.3057.1241703565.36',
            'modality': 'CT',
            'time': '000000'
        }
        self.assertEqual(self.dp.GetSeriesInfo(), series)

    def test_frame_of_reference(self):
        """Test if the frame of reference can be parsed."""
        referenceframe = '2.16.840.1.113662.2.12.0.3057.1241703565.36'
        self.assertEqual(self.dp.GetFrameOfReferenceUID(), referenceframe)

    def test_demographics(self):
        """Test if the demographics info can be parsed."""
        patient = {
            'name': 'boost^breast',
            'given_name': 'breast',
            'middle_name': '',
            'family_name': 'boost',
            'id': '123456',
            'gender': 'O',
            'birth_date': None
        }
        self.assertEqual(self.dp.GetDemographics(), patient)


class TestImage(unittest.TestCase):
    """Unit tests for Image Modality."""

    def setUp(self):
        """Setup files for Image modality testing."""
        ct_0_dcm = os.path.join(example_data, "ct.0.dcm")
        self.dp = dicomparser.DicomParser(ct_0_dcm)

    def test_image_data(self):
        """Test if the image data info can be parsed."""
        data = {
            'position': mv(DSfloat, ['-275', '-524', '168.5593']),
            'orientation':
                mv(DSfloat, ['1', '0.0', '-1.224647e-16', '0.0', '1', '0.0']),
            'pixelspacing': mv(DSfloat, ['1.074219', '1.074219']),
            'rows': 512,
            'columns': 512,
            'samplesperpixel': 1,
            'photometricinterpretation': 'MONOCHROME2',
            'littlendian': True,
            'patientposition': 'HFS',
            'frames': 1
        }
        self.assertEqual(self.dp.GetImageData(), data)

    def test_image_location(self):
        """Test if the image location can be parsed."""
        loc = 168.55929999999998
        self.assertEqual(self.dp.GetImageLocation(), loc)

    def test_image_orientation(self):
        """Test if the image orientation can be parsed."""
        orientation = 'SA'
        self.assertEqual(self.dp.GetImageOrientationType(), orientation)

    @unittest.skipUnless(pil_available, "PIL not installed")
    def test_image_generation(self):
        """Test if the image can be generated."""
        image = 90
        self.assertEqual(self.dp.GetImage().getpixel((255, 254)), image)

    def test_patient_to_pixel_lut(self):
        """Test if the image transformation matrix (LUT) can be generated."""
        lutvalue = 273.925909
        doselut = self.dp.GetPatientToPixelLUT()
        self.assertAlmostEqual(doselut[0][-1], lutvalue)


class TestRTStructureSet(unittest.TestCase):
    """Unit tests for RT Structure Set Modality."""

    def setUp(self):
        """Setup the files for RT Structure Set modality testing."""
        rtss_dcm = os.path.join(example_data, "rtss.dcm")
        self.dp = dicomparser.DicomParser(rtss_dcm)

    def test_referenced_series(self):
        """Test if the referenced series can be parsed."""
        series = '2.16.840.1.113662.2.12.0.3057.1241703565.43'
        self.assertEqual(self.dp.GetReferencedSeries(), series)

    def test_referenced_frame_of_reference(self):
        """Test if the referenced frame of reference can be parsed."""
        referenceframe = '2.16.840.1.113662.2.12.0.3057.1241703565.36'
        self.assertEqual(self.dp.GetFrameOfReferenceUID(), referenceframe)

    def test_structure_set_info(self):
        """Test if the structure set info can be parsed."""
        structure = {
            'label': 'CT_1',
            'date': '19010101',
            'time': '000000',
            'numcontours': 10
        }
        self.assertEqual(self.dp.GetStructureInfo(), structure)

    def test_structure_set_data(self):
        """Test if the structure set data can be parsed."""
        structure = {
            'color': array([255, 128, 0]),
            'empty': False,
            'id': 5,
            'name': 'Heart',
            'type': 'ORGAN'
        }
        structures = self.dp.GetStructures()
        # Pop the color numpy array
        assert_array_equal(structures[5].pop('color'),
                           structure.pop('color'))
        self.assertEqual(structures[5], structure)

    def test_structure_set_coordinates(self):
        """Test if the structure set coordinates can be parsed."""
        plane = {
            'type': 'CLOSED_PLANAR',
            'num_points': 166,
            'data': array((2.69, -313.11, -47.44))
        }
        planes = self.dp.GetStructureCoordinates(5)
        # Pop the planes coordinate numpy array
        assert_array_equal(array(planes['-47.440'][0].pop('data')[0]),
                           plane.pop('data'))
        self.assertEqual(planes['-47.440'][0], plane)

    def test_structure_plane_thickness(self):
        """Test if a structure plane thickness can be parsed."""
        planes = self.dp.GetStructureCoordinates(5)
        thickness = 3
        self.assertAlmostEqual(
            self.dp.CalculatePlaneThickness(planes), thickness)

    @unittest.skipUnless(shapely_available, "shapely not installed")
    def test_structure_volume(self):
        """Test if a structure volume can be calculated."""
        coords = self.dp.GetStructureCoordinates(5)
        volume = 432.84104445
        self.assertAlmostEqual(
            self.dp.CalculateStructureVolume(coords, 3), volume)

    @unittest.skipUnless(shapely_available, "shapely not installed")
    def test_structure_volume_holes(self):
        """Test if a structure volume with holes can be calculated."""
        coords = self.dp.GetStructureCoordinates(6)
        volume = 1995.1847937
        self.assertAlmostEqual(
            self.dp.CalculateStructureVolume(coords, 3), volume)


class TestRTPlan(unittest.TestCase):
    """Unit tests for RT Plan Modality."""

    def setUp(self):
        """Setup the files for RT Plan modality testing."""
        rtplan_dcm = os.path.join(example_data, "rtplan.dcm")
        self.dp = dicomparser.DicomParser(rtplan_dcm)

    def test_referenced_structureset(self):
        """Test if the referenced structure set can be parsed."""
        rtss = '1.2.246.352.71.4.320687012.3190.20090511122144'
        self.assertEqual(self.dp.GetReferencedStructureSet(), rtss)

    def test_plan_data(self):
        """Test if the plan data can be parsed."""
        data = {
            'label': 'B1',
            'date': '19010101',
            'time': '000000',
            'name': 'Breast',
            'rxdose': 1400,
            'brachy': False,
        }

        # Test DoseRefererenceSequence
        plandata = self.dp.GetPlan()
        self.assertEqual(plandata, data)

        # Test FractionGroupSequence
        del self.dp.ds.DoseReferenceSequence
        plandata = self.dp.GetPlan()
        data['name'] = ''
        self.assertEqual(plandata, data)

    def test_plan_beams_in_fraction(self):
        """Test if beams for a given fraction data can be parsed."""
        self.maxDiff = None
        beamdata = {'description': '', 'dose': 350.0, 'name': '3 RAO'}
        beams = self.dp.GetReferencedBeamsInFraction()
        self.assertEqual(beams[1], beamdata)


class TestRTDose(unittest.TestCase):
    """Unit tests for RT Dose Modality."""

    def setUp(self):
        """Setup the files for RT Dose modality testing."""
        rtdose_dcm = os.path.join(example_data, "rtdose.dcm")
        self.dp = dicomparser.DicomParser(rtdose_dcm)

    def test_referenced_plan(self):
        """Test if the referenced plan can be parsed."""
        rtplan = '1.2.246.352.71.5.320687012.24189.20090603083342'
        self.assertEqual(self.dp.GetReferencedRTPlan(), rtplan)

    def test_has_dvhs(self):
        """Test if DVHs exist in the dataset."""
        self.assertTrue(self.dp.HasDVHs())

    def test_dvh_data(self):
        """Test if the DVH data can be parsed."""
        dvh = self.dp.GetDVHs()[5]
        dvh.rx_dose = 14
        assert_array_almost_equal(dvh.bins, arange(0, 3.12, 0.01))
        self.assertEqual(dvh.dvh_type, 'cumulative')
        self.assertEqual(dvh.dose_units, 'Gy')
        self.assertEqual(dvh.volume_units, 'cm3')
        self.assertEqual(dvh.volume, 437.4623175026471)
        self.assertEqual(dvh.relative_dose().max, 22.214285714285555)
        self.assertEqual(dvh.relative_dose().min, 0.14285714285714285)
        self.assertEqual(dvh.relative_dose().mean, 4.590916280337491)

    def test_dose_grid(self):
        """Test if the dose grid can be parsed."""
        mean = 393.83944697514585
        # Test for plane that is listed in the GFOV
        self.assertEqual(self.dp.GetDoseGrid(-122.4407).mean(), mean)
        # Test for plane that doesn't exist in the dose grid
        assert_array_equal(self.dp.GetDoseGrid(-10000), array([]))

    def test_isodose_points(self):
        """Test if isodose points can be generated from the dose grid."""
        points = [(106, 20), (108, 20), (110, 20)]
        self.assertEqual(self.dp.GetIsodosePoints()[0:3], points)

    def test_isodose_points_memmap(self):
        """Test if isodose points can be generated via a memmapped dose grid."""
        points = [(106, 20), (108, 20), (110, 20)]
        dp = dicomparser.DicomParser(
            os.path.join(example_data, "rtdose.dcm"),
            memmap_pixel_array=True)
        self.assertEqual(dp.GetIsodosePoints()[0:3], points)

    def test_dose_data(self):
        """Test if the dose data can be parsed."""
        data = {
            'position': mv(
                DSfloat, ['-228.6541915', '-419.2444776', '-122.4407']),
            'orientation':
                mv(DSfloat, ['1', '0.0', '0.0', '0.0', '1', '0.0']),
            'pixelspacing': mv(DSfloat, [2.5, 2.5]),
            'rows': 129,
            'columns': 194,
            'samplesperpixel': 1,
            'photometricinterpretation': 'MONOCHROME2',
            'littlendian': True,
            'frames': 98,
            'doseunits': 'GY',
            'dosetype': 'PHYSICAL',
            'dosecomment': '',
            'dosesummationtype': 'PLAN',
            'dosegridscaling': 1.4e-05,
            'dosemax': 1048626.0,
            'lut': 253.8458085,
            'fraction': '',
            'x_lut_index': 0,
        }
        dosedata = self.dp.GetDoseData()
        # Pop the LUT numpy array
        assert_array_almost_equal(
            dosedata.pop('lut')[0][-1], data.pop('lut'))
        self.assertEqual(dosedata, data)
        dp = dicomparser.DicomParser(
            os.path.join(example_data, "rtdose.dcm"),
            memmap_pixel_array=True)
        dosedata = dp.GetDoseData()
        data['lut'] = 253.8458085
        assert_array_almost_equal(
            dosedata.pop('lut')[0][-1], data.pop('lut'))
        self.assertEqual(dosedata, data)


if __name__ == '__main__':
    import sys
    sys.exit(unittest.main())
